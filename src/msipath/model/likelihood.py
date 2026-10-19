"""Single altered allele and fixed repeat-length shift, marginalized per locus."""
from functools import lru_cache
import numpy as np
from scipy.special import logsumexp, betaln
from .common import require, normalize

@lru_cache(maxsize=32)
def fraction_quadrature(order,maximum,a,b):
    require(order>=2 and 0<maximum<=1 and a>0 and b>0,'Invalid MSI quadrature/prior')
    nodes,weights=np.polynomial.legendre.leggauss(order)
    fractions=(nodes+1)*maximum/2
    logw=np.log(weights)+(a-1)*np.log(fractions)+(b-1)*np.log1p(-fractions)-betaln(a,b)
    return fractions,logw-logsumexp(logw)


def length_kernel(length,period,size,stutter,indel):
    require(0<=stutter<1 and 0<=indel<1 and stutter+indel<1,'Invalid error probabilities')
    p=np.zeros(size); p[length]+=1-stutter-indel
    for shift in range(1,7):
        mass=stutter*(.5**shift)/(1-.5**6)/2
        for sign in (-1,1): p[min(size-1,max(0,length+sign*period*shift))]+=mass
    for sign in (-1,1): p[min(size-1,max(0,length+sign))]+=indel/2
    return normalize(p)


def normal_profile_key(locus):
    """Reference-only motif/length grouping, fixed before seeing query reads."""
    median=float(np.median([len(a['repeat_sequence']) for a in locus['alleles']]))
    bucket='short' if median<30 else 'medium' if median<60 else 'long'
    return locus['motif'].upper()+':'+bucket


def locus_settings(locus, settings):
    result=dict(settings)
    profiles=settings.get('normal_profiles',{})
    selected=profiles.get(normal_profile_key(locus),profiles.get('__global__',{}))
    allowed={'stutter_rate','indel_rate','novel_max_shift_repeats','novel_geometric_decay','out_of_panel_prior'}
    require(set(selected)<=allowed,'Invalid frozen normal profile keys')
    result.update(selected)
    return result


def structured_bf(locus, paths, observations, settings):
    settings = locus_settings(locus, settings)
    period = len(locus['motif'])
    lengths = [len(a['repeat_sequence']) for a in locus['alleles']]
    size = max(lengths) + period * 46 + 2
    if any(o['length_bp'] >= size for o in observations):
        raise ValueError('Observed tract outside fixed support')
    # Combine identical emissions; retain mapping weights exactly.
    keys = sorted({(int(o['length_bp']), float(o['map_weight'])) for o in observations})
    counts = np.array([sum((o['length_bp'], float(o['map_weight'])) == k
                           for o in observations) for k in keys])
    x = np.array([k[0] for k in keys]); weights = np.array([k[1] for k in keys])
    cache = {}

    def kernel(n):
        if n not in cache:
            p = length_kernel(n, period, size, settings['stutter_rate'], settings['indel_rate'])
            cache[n] = weights*p[x] + (1-weights)/size
        return cache[n]

    balances = settings['allele_balance_grid']
    bweights = np.asarray(settings['allele_balance_weights'], float)
    bweights /= bweights.sum()
    d = np.arange(1, 41)
    short = .5**(d-1); short /= short.sum()
    broad = .9**(d-1); broad /= broad.sum()
    shift_weights = .7*short + .3*broad
    normal_terms = []; directional = {-1: [], 1: []}
    gtotal = sum(g['weight'] for g in paths['genotypes'])
    for g in paths['genotypes']:
        if g['weight'] <= 0:
            continue
        a, b = lengths[g['allele1']], lengths[g['allele2']]
        ka, kb = kernel(a), kernel(b)
        for balance, bw in zip(balances, bweights):
            lp = np.log(g['weight']/gtotal) + np.log(bw)
            e0 = balance*ka + (1-balance)*kb
            normal_terms.append(lp + np.log(np.maximum(e0, 1e-300)) @ counts)
            for n, old, mass in ((a, ka, balance), (b, kb, 1-balance)):
                fractions, fw = fraction_quadrature(12, min(.5, mass), 1., 9.)
                for sign in (-1, 1):
                    allowed = [(int(v), float(w)) for v,w in zip(d, shift_weights)
                               if n+sign*int(v)*period >= period]
                    if not allowed:
                        directional[sign].append(lp-np.log(2)+np.log(np.maximum(e0,1e-300))@counts)
                        continue
                    sw = np.array([w for _,w in allowed]); sw /= sw.sum()
                    changed = np.array([kernel(n+sign*v*period) for v,_ in allowed])
                    # Remove mass from the mutated allele only.
                    emissions = e0[None,None,:] + fractions[None,:,None]*(changed[:,None,:]-old)
                    terms = np.log(np.maximum(emissions, 1e-300)) @ counts
                    directional[sign].append(lp-np.log(2)+logsumexp(terms+fw[None,:]+np.log(sw)[:,None]))
    l0 = logsumexp(normal_terms)
    contraction = logsumexp(directional[-1])-l0
    expansion = logsumexp(directional[1])-l0
    return dict(log_bf=float(np.logaddexp(contraction, expansion)-np.log(2)),
                contraction_log_bf=float(contraction), expansion_log_bf=float(expansion))


def aggregate(loci):
    """Scores only: dependence is handled by independent sample calibration."""
    if len(loci) < 5:
        return None
    b = np.array([r['log_bf'] for r in loci])
    # Fixed robust positive-evidence sum, identical for original/A/B.
    score = float(np.clip(b, 0, 10).sum())
    return score

