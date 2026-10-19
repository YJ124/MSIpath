"""Streaming SAM/BAM evidence with CIGAR-aware tracts and fragment accounting."""
from __future__ import annotations
import bisect
import re
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from .common import require, write_table, digest

def alignment_records(path, samtools='samtools', reference=None):
    path = Path(path)
    if path.suffix.lower()=='.sam':
        with path.open() as f:
            yield from f
    else:
        cmd = [samtools,'view','-h']
        if path.suffix.lower()=='.cram':
            require(reference,'CRAM requires the original reference FASTA')
            cmd += ['-T',str(reference)]
        cmd += [str(path)]
        with subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True) as p:
            yield from p.stdout
            err = p.stderr.read()
            require(p.wait()==0, f'samtools failed: {err}')

def cigar_map(start, cigar):
    """Map reference positions (0-based) to query indices; deletions are None."""
    tokens = re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    require(''.join(n+op for n,op in tokens)==cigar, f'Invalid CIGAR: {cigar}')
    r,q = start,0
    mapping = {}
    for count,op in tokens:
        n = int(count)
        if op in 'M=X':
            mapping.update((r+i,q+i) for i in range(n)); r+=n; q+=n
        elif op in 'DN':
            mapping.update((r+i,None) for i in range(n)); r+=n
        elif op in 'IS': q+=n
    return mapping,r,q

def tract_observation(seq, quality, mapping, locus, mapq, settings):
    start,end = locus['start']-1,locus['end']
    anchor = int(settings.get('anchor_bases',5))
    left = [mapping.get(p) for p in range(start-anchor,start)]
    right = [mapping.get(p) for p in range(end,end+anchor)]
    left_ok = all(q is not None for q in left)
    right_ok = all(q is not None for q in right)
    minq = float(settings.get('min_anchor_mean_quality',20))
    if left_ok: left_ok = np.mean([quality[q] for q in left])>=minq
    if right_ok: right_ok = np.mean([quality[q] for q in right])>=minq
    period = len(locus['motif'])
    if left_ok and right_ok:
        a,b = left[-1]+1,right[0]
        if b<a: return None
        return dict(kind='exact',length_bp=b-a,map_weight=1-10**(-mapq/10),mean_quality=float(np.mean(quality[a:b])) if b>a else minq)
    # One-boundary observations are right-censored length evidence, never
    # converted to complete lengths. Do not use unanchored in-repeat reads.
    if settings.get('use_censored_reads',True) and (left_ok or right_ok):
        positions = [q for p,q in mapping.items() if start<=p<end and q is not None]
        if not positions: return None
        a,b = (left[-1]+1,max(positions)+1) if left_ok else (min(positions),right[0])
        if b-a<2*period: return None
        return dict(kind='lower_bound',length_bp=b-a,map_weight=1-10**(-mapq/10),mean_quality=float(np.mean(quality[a:b])))
    return None

def collect_evidence(path, panel, settings, samtools='samtools', reference=None):
    by_chr = defaultdict(list)
    for i,s in enumerate(panel['variants']): by_chr[s['chrom']].append((s['pos']-1,i))
    snp_positions = {c:[x[0] for x in a] for c,a in by_chr.items()}
    loci = defaultdict(list)
    for loc,r in panel['loci'].items(): loci[r['chrom']].append((r['start']-1,r['end'],loc))
    for c in loci: loci[c].sort()
    locus_starts = {c:[x[0] for x in a] for c,a in loci.items()}
    max_length = {c:max(e-s for s,e,_ in a) for c,a in loci.items()}
    bases = defaultdict(dict); repeats = defaultdict(dict); qc = Counter()
    for line in alignment_records(path,samtools,reference):
        if line.startswith('@'): continue
        a = line.rstrip().split('\t')
        if len(a)<11: continue
        qc['alignment_records']+=1
        flag,mapq = int(a[1]),int(a[4])
        if flag & (4|256|512|1024|2048): qc['flag_excluded']+=1; continue
        if mapq<int(settings.get('min_mapq',10)): qc['mapq_excluded']+=1; continue
        if 'N' in a[5]: qc['reference_skip_excluded']+=1; continue
        chrom,start,seq = a[2],int(a[3])-1,a[9]
        if seq=='*' or a[10]=='*': qc['missing_sequence_quality']+=1; continue
        qualities = [ord(x)-33 for x in a[10]]
        require(len(seq)==len(qualities),'SAM sequence/quality length mismatch')
        mapping,end,query_len = cigar_map(start,a[5])
        require(query_len==len(seq),f'CIGAR/query length mismatch for {a[0]}')
        name=a[0]
        # Use read-group ID to avoid collapsing unrelated fragments with the
        # same name across sequencing libraries.
        group=next((x[5:] for x in a[11:] if x.startswith('RG:Z:')),'')
        fragment=(group,name)
        positions = snp_positions.get(chrom,[])
        for _,i in by_chr.get(chrom,[])[bisect.bisect_left(positions,start):bisect.bisect_left(positions,end)]:
            s=panel['variants'][i]; q=mapping.get(s['pos']-1)
            if q is None or qualities[q]<int(settings.get('min_snp_baseq',20)): continue
            base=seq[q].upper()
            if base not in (s['ref'],s['alt']): qc['snp_other_base']+=1; continue
            value=(int(base==s['alt']),qualities[q])
            old=bases[i].get(fragment)
            if old is None or value[1]>old[1]: bases[i][fragment]=value
            elif value[1]==old[1] and value[0]!=old[0]: bases[i][fragment]=(-1,value[1])
        starts=locus_starts.get(chrom,[])
        lo=bisect.bisect_left(starts,start-max_length.get(chrom,0))
        hi=bisect.bisect_left(starts,end)
        for ls,le,loc in loci.get(chrom,[])[lo:hi]:
            if le<=start: continue
            qc['locus_read_overlaps']+=1
            obs=tract_observation(seq,qualities,mapping,panel['loci'][loc],mapq,settings)
            if obs is None: qc['no_boundary_evidence']+=1; continue
            old=repeats[loc].get(fragment)
            if old is None: repeats[loc][fragment]=obs
            elif old.get('conflict'): continue
            elif old['kind']=='exact' and obs['kind']=='exact':
                if old['length_bp']!=obs['length_bp']:
                    repeats[loc][fragment]={'conflict':True}; qc['mate_length_conflicts']+=1
                else: old['map_weight']=max(old['map_weight'],obs['map_weight'])
            elif old['kind']=='lower_bound' and obs['kind']=='exact':
                if obs['length_bp']<old['length_bp']: repeats[loc][fragment]={'conflict':True}
                else: repeats[loc][fragment]=obs
            elif old['kind']=='exact' and obs['length_bp']>old['length_bp']:
                repeats[loc][fragment]={'conflict':True}
            elif old['kind']=='lower_bound' and obs['length_bp']>old['length_bp']: repeats[loc][fragment]=obs
    snp_counts={i:(sum(v[0]==0 for v in data.values()),sum(v[0]==1 for v in data.values())) for i,data in bases.items()}
    observations={loc:[dict(v,read_id=name,read_group=group) for (group,name),v in values.items() if not v.get('conflict')] for loc,values in repeats.items()}
    qc['snp_sites_observed']=sum(sum(x)>0 for x in snp_counts.values())
    qc['exact_fragments']=sum(o['kind']=='exact' for a in observations.values() for o in a)
    qc['censored_fragments']=sum(o['kind']=='lower_bound' for a in observations.values() for o in a)
    return snp_counts,observations,dict(qc)

def map_fastq(sample, directory, tools, fasta):
    """External aligner is a candidate-placement stage, followed by path inference."""
    d=Path(directory); d.mkdir(parents=True,exist_ok=True)
    r1,r2=Path(sample['r1']),Path(sample['r2'])
    require(r1.is_file() and r2.is_file(),'Missing FASTQ pair')
    bam=d/'alignment.bam'
    require(not bam.exists(),'Refusing stale BAM reuse')
    threads=int(tools.get('threads',4))
    bwa=[tools.get('bwa','bwa-mem2'),'mem','-t',str(threads),str(fasta),str(r1),str(r2)]
    sort=[tools.get('samtools','samtools'),'sort','-@',str(max(1,threads//2)),'-o',str(bam),'-']
    with (d/'mapping.log').open('w') as log:
        log.write(repr(bwa)+'\n'+repr(sort)+'\n'); log.flush()
        with subprocess.Popen(bwa,stdout=subprocess.PIPE,stderr=log) as p:
            result=subprocess.run(sort,stdin=p.stdout,stderr=log)
            p.stdout.close(); code=p.wait()
            require(code==0 and result.returncode==0,'FASTQ mapping failed; see mapping.log')
        subprocess.run([tools.get('samtools','samtools'),'index',str(bam)],stderr=log,check=True)
        subprocess.run([tools.get('samtools','samtools'),'quickcheck',str(bam)],stderr=log,check=True)
    return bam,dict(r1_sha256=digest(r1),r2_sha256=digest(r2),mapping_command=bwa)
