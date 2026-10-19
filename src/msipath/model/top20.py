"""SNP-only personalized reference selection and normal rank calibration."""
import copy
import numpy as np
from .common import require
from .haplotypes import ancestry, genotype_log_likelihood

def subset_panel(panel, donor_ids):
    keep=set(donor_ids)
    indices=[i for i,d in enumerate(panel['donors']) if d in keep]
    require(set(panel['donors'][i] for i in indices)==keep,'Missing requested donor')
    require(len(indices)==2*len(keep),'Expected two phased haplotypes per donor')
    p={k:v for k,v in panel.items() if k not in ('snps','str_haps','loci','phase_compatible')}
    for key in ('donors','hap_ids','populations'): p[key]=[panel[key][i] for i in indices]
    p['snps']=panel['snps'][:,indices].copy(); p['str_haps']={}; p['loci']={}; p['phase_compatible']={}
    populations=np.array(p['populations'])
    for locus in panel['locus_order']:
        old=panel['str_haps'][locus][indices]; used=sorted(set(int(v) for v in old if v>=0))
        lookup={a:b for b,a in enumerate(used)}
        h=np.array([lookup.get(int(v),-1) for v in old],dtype=np.int32)
        loc=copy.deepcopy(panel['loci'][locus]); loc['alleles']=[]
        for a in used:
            allele=copy.deepcopy(panel['loci'][locus]['alleles'][a])
            allele['allele_id']=f'a{lookup[a]}'
            allele['frequencies']={pop:float(np.sum((old==a)&(populations==pop))/max(1,np.sum((old>=0)&(populations==pop)))) for pop in set(populations)}
            loc['alleles'].append(allele)
        p['loci'][locus]=loc; p['str_haps'][locus]=h
        for key,allowed in panel['phase_compatible'].get(locus,{}).items():
            mapped='ALL' if allowed=='ALL' else [j for j,i in enumerate(indices) if i in set(allowed)]
            p['phase_compatible'].setdefault(locus,{})[key]=mapped
        p['phase_compatible'].setdefault(locus,{})
    return p


def recommend(panel,counts,settings,exclude_donor=None):
    donors=sorted(set(panel['donors'])-{exclude_donor})
    k=int(settings.get('top_n',20)); require(len(donors)>=k,'Insufficient reference donors')
    eligible=subset_panel(panel,donors)
    mixture=ancestry(eligible,counts,settings.get('ancestry',{'min_snps':30}))
    # Candidate ranking has equal contributions from mean normalized SNP
    # log-likelihood in chromosome-wide markers and local flanking markers.
    # Read-backed SNP phasing is not inferred here; phased copying is Step4.
    observed=[i for i in range(len(panel['variants'])) if sum(counts.get(i,(0,0)))>=4]
    global_ids=[i for i in observed if panel['variants'][i]['ancestry']]
    local_set={int(i) for values in panel['phase_compatible'].values() for i in values}
    local_ids=[i for i in observed if i in local_set]
    require(len(global_ids)>=30 and len(local_ids)>=3,'Insufficient recommendation SNPs')
    gl={}
    for i in set(global_ids+local_ids):
        value=genotype_log_likelihood(*counts[i]); gl[i]=(value-value.max())/sum(counts[i])
    rows=[]
    for d in donors:
        ix=[i for i,v in enumerate(eligible['donors']) if v==d]
        def score(ids):
            vals=[float(gl[i][int(eligible['snps'][i,ix].sum())]) for i in ids if np.all(eligible['snps'][i,ix]>=0)]
            require(len(vals)>=.8*len(ids),'Reference donor has excessive missing SNPs: '+d)
            return float(np.mean(vals)),len(vals)
        global_score,ng=score(global_ids); flank_score,nl=score(local_ids)
        pop=eligible['populations'][ix[0]]
        anc=float(np.log(max(mixture['proportions'].get(pop,0),.01)))
        total=global_score+flank_score+float(settings.get('ancestry_weight',.1))*anc
        rows.append(dict(donor=d,population=pop,score=total,global_snp_score=global_score,
                         flank_snp_score=flank_score,ancestry_log_compatibility=anc,
                         global_snps=ng,flank_snps=nl))
    rows.sort(key=lambda r:(-r['score'],r['donor']))
    for i,r in enumerate(rows,1): r.update(rank=i,selected=i<=k)
    selected=[r['donor'] for r in rows[:k]]
    require(exclude_donor not in selected,'Self-reference leakage')
    return selected,rows,mixture


def rank_p(score,null):
    require(len(null)>0 and np.all(np.isfinite(null)),'Invalid independent calibration scores')
    return (1+sum(v>=score-1e-12 for v in null))/(len(null)+1)

