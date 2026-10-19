"""Supervised admixture and diploid Li–Stephens-style copying HMM.

The HMM is exact within its explicitly reported donor subset. It is not a
wrapper around ADMIXTURE and does not claim to reproduce its optimizer.
"""
import math
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, xlogy
from .common import require, normalize

def genotype_log_likelihood(ref_count,alt_count,error=.005):
    f=np.array([error,.5,1-error])
    return alt_count*np.log(f)+ref_count*np.log1p(-f)

def ancestry(panel, counts, settings):
    labels=sorted(set(panel['populations']))
    pops=np.array(panel['populations'])
    indices=[i for i,s in enumerate(panel['variants']) if s['ancestry'] and sum(counts.get(i,(0,0)))>0]
    require(len(indices)>=int(settings.get('min_snps',30)),f'INSUFFICIENT_ANCESTRY_SNPS: {len(indices)} observed')
    freq=[]; gl=[]; valid=[]
    for i in indices:
        row=panel['snps'][i]
        f=[]
        for pop in labels:
            a=row[pops==pop]; a=a[a>=0]
            f.append((a.sum()+.5)/(len(a)+1))
        freq.append(f); gl.append(genotype_log_likelihood(*counts[i])); valid.append(i)
    freq=np.array(freq); gl=np.array(gl)
    def objective(q):
        f=np.clip(freq@q,1e-6,1-1e-6)
        prior=np.column_stack(((1-f)**2,2*f*(1-f),f*f))
        return -float(logsumexp(gl+np.log(prior),axis=1).sum())
    if len(labels)==1: q=np.ones(1); success=True
    else:
        result=minimize(objective,np.full(len(labels),1/len(labels)),method='SLSQP',bounds=[(0,1)]*len(labels),
                        constraints=[{'type':'eq','fun':lambda x:x.sum()-1}],options={'maxiter':500,'ftol':1e-9})
        q=normalize(np.maximum(result.x,0)); success=result.success
    require(success,'Ancestry mixture optimizer did not converge')
    # PCA projection is descriptive; genotype likelihoods, not PCA coordinates,
    # determine the mixture used for the copying prior.
    h=panel['snps'][valid].astype(float)
    means=np.array([np.mean(a[a>=0]) if np.any(a>=0) else .5 for a in h])
    h=np.where(h<0,means[:,None],h)
    dosage=h[:,0::2]+h[:,1::2]
    center=dosage.mean(axis=1); scale=np.sqrt(np.maximum(2*means*(1-means),.02))
    x=((dosage-center[:,None])/scale[:,None]).T
    _,_,vh=np.linalg.svd(x,full_matrices=False)
    f=np.clip(freq@q,1e-6,1-1e-6)
    post=gl+np.log(np.column_stack(((1-f)**2,2*f*(1-f),f*f)))
    post=np.exp(post-logsumexp(post,axis=1)[:,None])
    target=post@np.arange(3)
    pcs=vh[:2]@((target-center)/scale)
    return dict(proportions=dict(zip(labels,map(float,q))),method='supervised_genotype_likelihood_mixture',
                snps=len(valid),observed_chromosomes=sorted({panel['variants'][i]['chrom'] for i in valid}),pca=pcs.tolist())

def forward_transition(a,w,rho):
    """T^T A T with independent copying transitions, O(H²)."""
    stay=1-rho
    return stay*stay*a+rho*stay*(w[:,None]*a.sum(axis=0)[None,:]+a.sum(axis=1)[:,None]*w[None,:])+rho*rho*np.outer(w,w)*a.sum()

def backward_transition(b,w,rho):
    stay=1-rho
    return stay*stay*b+rho*stay*((w@b)[None,:]+(b@w)[:,None])+rho*rho*(w@b@w)

def pair_emission(h,counts,epsilon):
    # Missing donor SNPs have a 0.5 allele probability; no arbitrary phase call.
    h=np.where(h<0,.5,h)
    p0=np.outer(1-h,1-h); p2=np.outer(h,h); p1=1-p0-p2
    g=genotype_log_likelihood(*counts,error=epsilon); g-=g.max()
    return np.maximum(p0*np.exp(g[0])+p1*np.exp(g[1])+p2*np.exp(g[2]),1e-300)

def match_paths(panel,locus,counts,mixture,settings):
    loc=panel['loci'][locus]; alleles=panel['str_haps'][locus]
    valid=np.flatnonzero(alleles>=0)
    require(len(valid)>=2,'NO_PHASED_STR_DONORS')
    pop=np.array(panel['populations']); weights=np.array([mixture.get(p,0)/max(1,sum(pop==p)) for p in pop])
    # WebSTR/registry population frequencies reweight the observed donor panel
    # by frequency-ratio, avoiding a second multiplication of the same prior.
    for p in set(pop):
        donor_mask=(pop==p)&(alleles>=0)
        n=donor_mask.sum()
        for a in np.unique(alleles[donor_mask]):
            target=loc['alleles'][int(a)].get('frequencies',{}).get(p)
            if target is not None:
                observed=sum(donor_mask&(alleles==a))/max(n,1)
                weights[donor_mask&(alleles==a)]*=max(float(target),1e-8)/observed
    available=[]
    compatible=panel['phase_compatible'].get(locus,{})
    for key,allowed in compatible.items():
        i=int(key)
        candidates=valid if allowed=='ALL' else np.intersect1d(valid,np.asarray(allowed,dtype=int))
        observed=panel['snps'][i,candidates]
        if sum(counts.get(i,(0,0)))>0 and len(candidates)>=2 and len(np.unique(observed[observed>=0]))>1: available.append(i)
    available.sort(key=lambda i:panel['variants'][i]['pos'])
    require(len(available)>=int(settings.get('min_anchor_snps',3)),f'INSUFFICIENT_LOCAL_SNPS:{len(available)}')
    # Preselect using SNP evidence ONLY; STR repeat reads must not influence the
    # normal-path prior. The cap and retained IDs are written into the outputs.
    scores=np.log(np.maximum(weights[valid],1e-300))
    epsilon=float(settings.get('snp_error',.005))
    for i in available:
        hs=panel['snps'][i,valid].astype(float)
        allowed=compatible[str(i)]
        if allowed!='ALL': hs[~np.isin(valid,allowed)]=-1
        hs=np.where(hs<0,.5,hs)
        g=genotype_log_likelihood(*counts[i],epsilon)
        # Integrate the unknown other chromosome using population allele mean.
        af=np.average(hs,weights=np.maximum(weights[valid],1e-300))
        emis=(1-hs)*((1-af)*np.exp(g[0]-g.max())+af*np.exp(g[1]-g.max()))+hs*((1-af)*np.exp(g[1]-g.max())+af*np.exp(g[2]-g.max()))
        scores+=np.log(np.maximum(emis,1e-300))
    cap=int(settings.get('max_haplotypes',128))
    require(cap>=2,'At least two donor candidates are required')
    ranked=np.argsort(-scores,kind='stable')
    if len(valid)<=cap:
        selected=valid[ranked]
    else:
        # Keep complementary chromosomes: selecting only the individually best
        # donors can drop one side of an otherwise obvious heterozygous pair.
        complementary=np.full(len(valid),-np.inf)
        for seed in ranked[:min(4,len(ranked))]:
            pair_score=np.log(np.maximum(weights[valid],1e-300))
            for i in available:
                h=panel['snps'][i,valid].astype(float)
                allowed=compatible[str(i)]
                if allowed!='ALL': h[~np.isin(valid,allowed)]=-1
                h=np.where(h<0,.5,h); seed_allele=h[seed]
                p0=(1-h)*(1-seed_allele); p2=h*seed_allele; p1=1-p0-p2
                g=genotype_log_likelihood(*counts[i],epsilon); g-=g.max()
                pair_score+=np.log(np.maximum(p0*np.exp(g[0])+p1*np.exp(g[1])+p2*np.exp(g[2]),1e-300))
            complementary=np.maximum(complementary,pair_score)
        keep=list(ranked[:cap//2])
        for j in np.argsort(-complementary,kind='stable'):
            if j not in keep: keep.append(int(j))
            if len(keep)==cap: break
        selected=valid[keep]
    w=normalize(weights[selected])
    emissions=[]; positions=[]
    center=(loc['start']+loc['end'])/2
    events=[(panel['variants'][i]['pos'],i) for i in available]+[(center,None)]
    events.sort(key=lambda x:x[0])
    for pos,i in events:
        if i is None: em=np.ones((len(selected),len(selected)))
        else:
            h=panel['snps'][i,selected].astype(float)
            allowed=compatible[str(i)]
            if allowed!='ALL': h[~np.isin(selected,allowed)]=-1
            em=pair_emission(h,counts[i],epsilon)
        emissions.append(em); positions.append(pos)
    target=next(i for i,(_,v) in enumerate(events) if v is None)
    rate=float(settings.get('switch_rate_per_bp',1e-6))
    rhos=[0]+[1-math.exp(-rate*(positions[i]-positions[i-1])) for i in range(1,len(positions))]
    a=np.outer(w,w); forward=[]
    for em,rho in zip(emissions,rhos):
        a=forward_transition(a,w,rho)*em; a/=a.sum(); forward.append(a)
    b=np.ones_like(a)
    for i in range(len(emissions)-1,target,-1):
        b=backward_transition(emissions[i]*b,w,rhos[i]); b/=b.max()
    posterior=forward[target]*b; posterior/=posterior.sum()
    pairs=[]; genotype_weights={}
    for i in range(len(selected)):
        for j in range(i,len(selected)):
            weight=float(posterior[i,j]+(posterior[j,i] if i!=j else 0))
            ai,aj=int(alleles[selected[i]]),int(alleles[selected[j]])
            key=tuple(sorted((ai,aj)))
            genotype_weights[key]=genotype_weights.get(key,0)+weight
            pairs.append(dict(hap1=panel['hap_ids'][selected[i]],hap2=panel['hap_ids'][selected[j]],allele1=ai,allele2=aj,weight=weight))
    pairs.sort(key=lambda x:-x['weight'])
    top=pairs[:int(settings.get('top_k',20))]
    return dict(topk=top,topk_mass=sum(r['weight'] for r in top),genotypes=[dict(allele1=a,allele2=b,weight=v) for (a,b),v in genotype_weights.items()],
        selected_haplotypes=selected.tolist(),eligible_haplotypes=len(valid),anchor_snps=len(available),
        donor_subset_approximation=len(selected)<len(valid))

def export_local_graph(path,panel,locus,paths):
    """Emit phased donor paths through SNP-bearing flanks and STR sequence nodes."""
    loc=panel['loci'][locus]; hap_lookup={h:i for i,h in enumerate(panel['hap_ids'])}
    used=sorted({p['hap1'] for p in paths['topk']}|{p['hap2'] for p in paths['topk']})
    nodes={}; links=set(); output_paths=[]
    for hap in used:
        h=hap_lookup[hap]; ai=int(panel['str_haps'][locus][h])
        left=list(loc['left_flank']); right=list(loc['right_flank'])
        left_start=loc['start']-len(left)
        for i,s in enumerate(panel['variants']):
            if s['chrom']!=loc['chrom']: continue
            allele=panel['snps'][i,h]
            allowed=panel['phase_compatible'].get(locus,{}).get(str(i),[])
            if allele<0 or (allowed!='ALL' and h not in allowed): continue
            base=s['alt'] if allele==1 else s['ref']
            if left_start<=s['pos']<loc['start']: left[s['pos']-left_start]=base
            if loc['end']<s['pos']<=loc['end']+len(right): right[s['pos']-loc['end']-1]=base
        ids=[]
        for typ,seq in [('L',''.join(left)),('A',loc['alleles'][ai]['repeat_sequence']),('R',''.join(right))]:
            key=(typ,seq)
            if key not in nodes: nodes[key]=f'n{len(nodes)+1}'
            ids.append(nodes[key])
        links.update(zip(ids,ids[1:])); output_paths.append((hap,ids))
    with open(path,'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for (typ,seq),nid in nodes.items(): f.write(f'S\t{nid}\t{seq or "*"}\tLN:i:{len(seq)}\tTP:Z:{typ}\n')
        for a,b in sorted(links): f.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
        for hap,ids in output_paths: f.write(f'P\t{hap}\t'+','.join(x+'+' for x in ids)+'\t*\n')
