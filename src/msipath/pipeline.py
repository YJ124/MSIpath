"""Default MSIpath workflow: fresh evidence, personalized paths, single-allele evidence."""
import hashlib
import json
import math
from pathlib import Path
import time
from . import __version__
from .model.common import ModelError, require, save_json, digest
from .model.reference import load_panel
from .model.evidence import collect_evidence
from .model.top20 import recommend, subset_panel, rank_p
from .model.haplotypes import match_paths, export_local_graph
from .model.likelihood import structured_bf, aggregate

DEFAULTS = {
    'evidence': {'min_mapq':10,'min_snp_baseq':20,'anchor_bases':5,'min_anchor_mean_quality':20,'use_censored_reads':True},
    'recommendation': {'top_n':20,'ancestry_weight':.1,'ancestry':{'min_snps':30}},
    'haplotypes': {'flank_bp':5000,'min_anchor_snps':3,'max_haplotypes':40,'top_k':20,'snp_error':.005,'switch_rate_per_bp':1e-6},
    'likelihood': {'stutter_rate':.03,'indel_rate':.01,'allele_balance_grid':[.25,.5,.75],
                   'allele_balance_weights':[.2,.6,.2],'max_length_bp':5000},
}

def canonical_hash(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def algorithm_hash():
    root=Path(__file__).parent
    files=[root/'pipeline.py']+sorted((root/'model').glob('*.py'))
    return canonical_hash({p.name:digest(p) for p in files})

def signature(panel_dir, profiles, loci):
    return dict(algorithm='single_allele_logBF_sum',version=__version__,code_sha256=algorithm_hash(),
        panel_sha256=digest(Path(panel_dir)/'panel.json'),profile_sha256=canonical_hash(profiles),
        settings_sha256=canonical_hash(DEFAULTS),loci_sha256=canonical_hash(loci))

def decide(score, bank, signature_value, donor, min_normals=20):
    if bank is not None:
        require(bank.get('signature')==signature_value,'Calibration does not match this model, panel, profiles or locus set')
        require(donor not in bank.get('donors',[]),'Query donor overlaps calibration donors')
    if score is None: return dict(call='UNCALLABLE',calibrated_p=None,calibration_n=0)
    require(math.isfinite(score),'Nonfinite sample score')
    if bank is None: return dict(call='UNCALIBRATED',calibrated_p=None,calibration_n=0)
    null=bank.get('scores',[])
    require(len(null)>=min_normals and all(math.isfinite(x) for x in null),'Insufficient or invalid calibration normals')
    p=rank_p(score,null)
    return dict(call='MSI' if p<=.05 else 'MSS',calibrated_p=p,calibration_n=len(null))

def infer(*, panel_dir, profile_file, output, sample_id, donor_id, bam=None,
          reference=None, evidence_dir=None, loci_file=None, calibration=None):
    start=time.monotonic(); out=Path(output)
    require(not out.exists(),'Output directory already exists; use a new directory')
    require(bool(bam)!=bool(evidence_dir),'Provide exactly one of BAM/CRAM/SAM or saved evidence')
    panel=load_panel(panel_dir); profiles=json.loads(Path(profile_file).read_text())
    wanted=json.loads(Path(loci_file).read_text()) if loci_file else list(panel['locus_order'])
    require(isinstance(wanted,list) and wanted and len(wanted)==len(set(wanted)),'Locus list must be nonempty and unique')
    require(set(wanted)<=set(panel['locus_order']),'Unknown locus in requested scope')
    loci=[l for l in panel['locus_order'] if l in set(wanted)]
    sig=signature(panel_dir,profiles,loci)
    bank=json.loads(Path(calibration).read_text()) if calibration else None
    # Reject calibration/scope errors before expensive BAM processing.
    decide(None,bank,sig,donor_id)
    out.mkdir(parents=True); (out/'graphs').mkdir()
    if bam:
        counts,obs,qc=collect_evidence(bam,panel,DEFAULTS['evidence'],reference=reference)
    else:
        evidence=Path(evidence_dir)
        counts={int(k):v for k,v in json.loads((evidence/'snp_counts.json').read_text()).items()}
        obs=json.loads((evidence/'repeat_fragments.json').read_text())
        qc={'input':'saved evidence; BAM extraction excluded'}
    save_json(out/'snp_counts.json',counts);save_json(out/'repeat_fragments.json',obs);save_json(out/'qc.json',qc)
    try:
        selected,ranking,mix=recommend(panel,counts,DEFAULTS['recommendation'],exclude_donor=donor_id)
    except ModelError as e:
        result=dict(method='MSIpath',version=__version__,sample_id=sample_id,donor=donor_id,
            score=None,callable_loci=0,signature=sig,call='UNCALLABLE',calibrated_p=None,calibration_n=0,
            reason=str(e),reference_donors=sorted(set(panel['donors'])),seconds=time.monotonic()-start)
        save_json(out/'result.json',result);return result
    local=subset_panel(panel,selected)
    save_json(out/'recommendation.json',dict(selected=selected,ranking=ranking,ancestry=mix))
    settings=dict(DEFAULTS['likelihood'],normal_profiles=profiles)
    results=[]; paths=[]
    for loc in loci:
        try:
            path=match_paths(local,loc,counts,mix['proportions'],DEFAULTS['haplotypes'])
            observations=obs.get(loc,[]); exact=[o for o in observations if o['kind']=='exact']
            require(len(exact)>=3 and len(observations)>=5,'INSUFFICIENT_FRAGMENTS')
            size=max(len(a['repeat_sequence']) for a in local['loci'][loc]['alleles'])+len(local['loci'][loc]['motif'])*46+2
            require(size<=settings['max_length_bp'] and all(o['length_bp']<size for o in observations),'OBSERVED_LENGTH_OUTSIDE_MODEL_SUPPORT')
            value=structured_bf(local['loci'][loc],path,exact,settings)
            require(all(math.isfinite(x) for x in value.values()),'NONFINITE_LIKELIHOOD')
            results.append(dict(locus=loc,status='PASS',exact_fragments=len(exact),**value))
            paths.append(dict(locus=loc,genotypes=path['genotypes'],anchor_snps=path['anchor_snps']))
            export_local_graph(out/'graphs'/(loc.replace(':','_')+'.gfa'),local,loc,path)
        except ModelError as e: results.append(dict(locus=loc,status=str(e)))
    passing=[r for r in results if r['status']=='PASS']; score=aggregate(passing)
    result=dict(method='MSIpath',version=__version__,sample_id=sample_id,donor=donor_id,score=score,
        callable_loci=len(passing),requested_loci=len(loci),signature=sig,reference_donors=sorted(set(panel['donors'])),
        input_mode='alignment' if bam else 'saved_evidence',**decide(score,bank,sig,donor_id))
    result['seconds']=time.monotonic()-start
    save_json(out/'loci.json',results);save_json(out/'paths.json',paths);save_json(out/'result.json',result)
    return result

def make_calibration(records, output):
    """records: explicitly declared independent normal donors and result.json paths."""
    require(not Path(output).exists(),'Calibration file already exists')
    require(len(records)>=20,'At least 20 independent normal donors are required')
    donors=[r['donor'] for r in records]
    require(len(set(donors))==len(donors),'Duplicate calibration donor')
    results=[]
    for r in records:
        require(r.get('role')=='calibration' and str(r.get('label'))=='0','Calibration requires declared normal calibration samples')
        value=json.loads(Path(r['result']).read_text())
        require(value['donor']==r['donor'],'Calibration donor/result mismatch')
        require(r['donor'] not in value['reference_donors'],'Calibration donor is in the reference panel')
        require(value['score'] is not None and math.isfinite(value['score']),'Unscorable calibration donor; supply sufficient usable independent normals')
        results.append(value)
    sig=results[0]['signature'];require(all(v['signature']==sig for v in results),'Calibration results use different model settings')
    bank=dict(method='MSIpath',signature=sig,donors=donors,scores=[v['score'] for v in results],
              normal_count=len(results),rule='upper-tail rank with plus one; alpha=0.05')
    Path(output).parent.mkdir(parents=True,exist_ok=True);save_json(output,bank);return bank

def tumor_purity_to_molecular_fraction(p):
    """Simulation conversion: clonal heterozygous alteration, diploid copy-neutral cells."""
    require(math.isfinite(p) and 0<=p<=1,'Tumor purity must be between zero and one')
    return p/2
