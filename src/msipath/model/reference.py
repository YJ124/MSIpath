"""Build a validated SNP/STR-linked panel from jointly phased sequence VCFs.

Coordinates: catalog start/end are 1-based inclusive; SNP pos is 1-based.
Independent SNP and STR phasing is deliberately NOT auto-joined by sample ID.
"""
from __future__ import annotations
import bisect
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
from .common import ModelError, require, text_open, read_table, save_json, digest, Fasta

def info_dict(text):
    return dict(x.split('=',1) if '=' in x else (x,True) for x in text.split(';') if x)

def build_panel(config, destination):
    destination = Path(destination)
    require(not destination.exists(), f'Panel already exists: {destination}')
    meta = read_table(config['sample_metadata'])
    donors = {r['sample_id']:r[config.get('population_column','super_population' if 'super_population' in r else 'population')] for r in meta}
    require(len(donors) == len(meta) and donors, 'Duplicated or empty donor metadata')
    excluded = set(config.get('exclude_samples', []))
    donors = {k:v for k,v in donors.items() if k not in excluded}
    phase = config.get('phase_provenance', {})
    require(phase.get('scope') == 'joint_snp_str_chromosome' and phase.get('source'),
            'Require documented jointly phased SNP–STR chromosome source; independent VCF GT order is not linkage')
    loci = {}
    pointers = set()
    inputs = list(config['gfa']) + list(config['registries']) + list(config['joint_phased_vcfs']) + [config['sample_metadata'],config['loci'],config['reference_fasta']]
    fasta=Fasta(config['reference_fasta'])
    for path in config['gfa']:
        with open(path) as f:
            for line in f:
                a = line.rstrip().split('\t')
                if a[:2] == ['X','PTR']: pointers.add(a[2])
                elif a[:2] == ['X','STR'] and config.get('allow_str_locus_pointers',False):
                    # Explicit locus-anchor adapter, not an assertion that every
                    # VCF sequence allele exists as a path in the source GFA.
                    pointers.add(a[2])
                elif a[0] == 'S': pointers.update(t[5:] for t in a[3:] if t.startswith('PI:Z:'))
    catalog = read_table(config['loci'])
    coord = {}
    for r in catalog:
        k = (r['chrom'],int(r['start']),int(r['end']))
        require(k not in coord, f'Duplicate locus coordinates: {k}')
        coord[k] = r['locus_id']
    for path in config['registries']:
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                key = r['chrom'],int(r['start']),int(r['end'])
                if key not in coord: continue
                require(r['pointer_id'] in pointers, f'Unknown GFA pointer {r["pointer_id"]}')
                require(coord[key] not in loci, f'Duplicate registry {key}')
                r['locus_id'] = coord[key]
                require(fasta.fetch(r['chrom'],r['start']-1,r['end'])==r['metadata']['reference_repeat_sequence'],f'FASTA/registry coordinate mismatch: {key}')
                r['alleles'] = [a for a in r['alleles'] if a.get('repeat_sequence') is not None]
                require(r['alleles'], f'No sequence alleles at {key}')
                loci[r['locus_id']] = r
    require(len(loci) == len(coord), 'Some requested loci lack graph/registry support')
    for loc,r in loci.items():
        for pop in set(donors.values()):
            require(any(float(a.get('frequencies',{}).get(pop,0))>0 for a in r['alleles']),f'Missing population frequency prior: {loc}/{pop}')
    donor_names = sorted(donors)
    hap_ids = [f'{s}|{i}' for s in donor_names for i in (1,2)]
    populations = [donors[s] for s in donor_names for _ in (1,2)]
    flank = int(config.get('flank_bp',500000))
    windows = defaultdict(list)
    for r in loci.values(): windows[r['chrom']].append((r['start']-flank,r['end']+flank))
    ancestry = set()
    if config.get('ancestry_markers'):
        ancestry = {(r['chrom'],int(r['pos'])) for r in read_table(config['ancestry_markers'])}
        inputs.append(config['ancestry_markers'])
    require(ancestry, 'An independently LD-pruned genome-wide ancestry-marker list is required')
    variants, snp_rows, seen = [], [], set()
    str_haps = {loc: [-1]*len(hap_ids) for loc in loci}
    phase_blocks = defaultdict(dict)
    snp_blocks = []
    for path in config['joint_phased_vcfs']:
        with text_open(path) as f:
            header = None
            for line in f:
                if line.startswith('##'): continue
                if line.startswith('#CHROM'):
                    a = line.rstrip().split('\t')
                    header = {s:i for i,s in enumerate(a[9:])}
                    require(all(s in header for s in donor_names), 'Joint VCF missing panel donors')
                    continue
                if line.startswith('#') or not line.strip(): continue
                require(header is not None,'VCF header missing')
                a = line.rstrip().split('\t')
                chrom,pos,ref,alt = a[0],int(a[1]),a[3],a[4].split(',')
                if a[6] not in ('PASS','.'): continue
                info = info_dict(a[7])
                start = int(info.get('START',pos)); end = int(info.get('END',pos+len(ref)-1))
                loc = coord.get((chrom,start,end))
                is_snp = len(ref)==1 and len(alt)==1 and len(alt[0])==1 and set(ref+alt[0]) <= set('ACGT')
                if loc is None and not (is_snp and ((chrom,pos) in ancestry or any(x<=pos<=y for x,y in windows[chrom]))): continue
                fmt = a[8].split(':'); fields = [dict(zip(fmt,a[9+header[s]].split(':'))) for s in donor_names]
                genotypes, ps = [], []
                for field in fields:
                    gt = field.get('GT','.')
                    bits = gt.replace('|','/').split('/')
                    if len(bits)!=2 or '.' in bits or ('|' not in gt and bits[0]!=bits[1]):
                        genotypes.extend([-1,-1])
                    else: genotypes.extend(map(int,bits))
                    ps.extend([field.get('PS','CHROMOSOME')]*2)
                if loc is not None:
                    r = loci[loc]
                    left_trim = start-pos; right_trim = pos+len(ref)-1-end
                    require(left_trim>=0 and right_trim>=0, f'Invalid STR VCF padding at {loc}')
                    sequences = [ref]+alt
                    body = lambda s: s[left_trim:len(s)-right_trim if right_trim else len(s)]
                    lookup = {x['repeat_sequence']:i for i,x in enumerate(r['alleles'])}
                    require(body(ref) == r['metadata']['reference_repeat_sequence'], f'VCF/reference tract mismatch at {loc}')
                    require(all(set(s)<=set('ACGTN') for s in sequences), f'Symbolic STR allele at {loc}; normalize explicitly')
                    for h,g in enumerate(genotypes):
                        if g>=0 and g<len(sequences) and body(sequences[g]) in lookup:
                            v = lookup[body(sequences[g])]
                            require(str_haps[loc][h] in (-1,v),f'Conflicting STR phase at {loc}')
                            str_haps[loc][h] = v
                            phase_blocks[loc][h] = ps[h]
                elif is_snp:
                    require(fasta.fetch(chrom,pos-1,pos)==ref,f'SNP REF mismatch {chrom}:{pos}')
                    key = chrom,pos
                    require(key not in seen,f'Duplicate SNP position {key}')
                    seen.add(key)
                    variants.append(dict(chrom=chrom,pos=pos,ref=ref,alt=alt[0],ancestry=key in ancestry))
                    snp_rows.append(genotypes); snp_blocks.append(ps)
                    require(len(variants)<=int(config.get('max_panel_snps',200000)), 'Panel exceeds max_panel_snps; provide region-filtered VCFs, do not silently downsample')
    require(variants,'No SNPs found in the joint panel')
    order = sorted(range(len(variants)),key=lambda i:(variants[i]['chrom'],variants[i]['pos']))
    variants = [variants[i] for i in order]
    snps = np.array([snp_rows[i] for i in order],dtype=np.int8)
    blocks = [snp_blocks[i] for i in order]
    require(np.all((snps>=-1)&(snps<=1)), 'Non-biallelic SNP coding')
    # Preserve block compatibility per locus rather than pretending that | implies
    # phase across arbitrary PS boundaries.
    compatible = {}
    for loc,r in loci.items():
        require(sum(x>=0 for x in str_haps[loc])>=2, f'No jointly phased STR haplotypes at {loc}')
        indices = [i for i,s in enumerate(variants) if s['chrom']==r['chrom'] and r['start']-flank<=s['pos']<=r['end']+flank and not r['start']<=s['pos']<=r['end']]
        compatible[loc]={}
        eligible=set(phase_blocks[loc])
        for i in indices:
            allowed=[h for h in eligible if blocks[i][h]==phase_blocks[loc][h]]
            compatible[loc][str(i)]='ALL' if len(allowed)==len(eligible) else sorted(allowed)
    destination.mkdir(parents=True)
    np.savez_compressed(destination/'haplotypes.npz',snps=snps,**{f'str_{i}':np.array(str_haps[loc],dtype=np.int32) for i,loc in enumerate(loci)})
    panel = dict(schema_version=2,assembly='GRCh38',hap_ids=hap_ids,donors=[s for s in donor_names for _ in (1,2)],populations=populations,
        variants=variants,loci=loci,locus_order=list(loci),phase_compatible=compatible,phase_provenance=phase,excluded_samples=sorted(excluded),
        sources=[dict(path=str(Path(p).resolve()),sha256=digest(p)) for p in inputs],flank_bp=flank)
    panel['haplotype_data_sha256']=digest(destination/'haplotypes.npz')
    save_json(destination/'panel.json',panel)
    fasta.close()
    return panel

def load_panel(directory):
    p = Path(directory)
    panel = json.loads((p/'panel.json').read_text(encoding='utf-8'))
    require(panel.get('schema_version')==2,'Not a path-model reference panel')
    require(panel.get('haplotype_data_sha256')==digest(p/'haplotypes.npz'),'Reference haplotype checksum mismatch')
    with np.load(p/'haplotypes.npz',allow_pickle=False) as a:
        panel['snps'] = a['snps']
        panel['str_haps'] = {loc:a[f'str_{i}'] for i,loc in enumerate(panel['locus_order'])}
    require(panel['snps'].shape==(len(panel['variants']),len(panel['hap_ids'])),'SNP panel dimensions mismatch')
    for loc,h in panel['str_haps'].items():
        require(len(h)==len(panel['hap_ids']) and np.all((h>=-1)&(h<len(panel['loci'][loc]['alleles']))),'STR panel dimensions/codes invalid')
    return panel
