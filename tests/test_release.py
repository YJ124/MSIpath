import json,math
from pathlib import Path
import numpy as np
import pytest
from msipath.pipeline import decide,make_calibration,tumor_purity_to_molecular_fraction
from msipath.model.common import ModelError
from msipath.model.likelihood import structured_bf,aggregate
from msipath.model.top20 import subset_panel,recommend

def test_calibration_ties_and_scope_and_donor_leakage():
 sig={'algorithm':'single_allele_logBF_sum'}
 bank={'signature':sig,'scores':[2.]*40,'donors':['normal'+str(i) for i in range(40)]}
 assert decide(2,bank,sig,'test')['call']=='MSS'
 assert decide(3,bank,sig,'test')['calibrated_p']==1/41
 assert decide(3,None,sig,'test')['call']=='UNCALIBRATED'
 assert decide(None,bank,sig,'test')['call']=='UNCALLABLE'
 with pytest.raises(ModelError):decide(3,bank,{},'test')
 with pytest.raises(ModelError):decide(3,bank,sig,'normal0')

def test_single_allele_symmetry_and_positive_evidence():
 locus={'motif':'AC','alleles':[{'repeat_sequence':'AC'*10},{'repeat_sequence':'AC'*12}]}
 cfg={'stutter_rate':.03,'indel_rate':.001,'allele_balance_grid':[.25,.5,.75],'allele_balance_weights':[.2,.6,.2]}
 obs=[{'kind':'exact','length_bp':n,'map_weight':1.} for n in [20]*14+[24]*10+[28]*6]
 a=structured_bf(locus,{'genotypes':[{'allele1':0,'allele2':1,'weight':1.}]},obs,cfg)
 b=structured_bf(locus,{'genotypes':[{'allele1':1,'allele2':0,'weight':1.}]},obs,cfg)
 assert a['log_bf']==pytest.approx(b['log_bf'],abs=1e-10)
 assert aggregate([a]*4) is None
 assert aggregate([{'log_bf':-2},{'log_bf':0},{'log_bf':2},{'log_bf':20},{'log_bf':4}])==16

def test_purity_conversion():
 assert tumor_purity_to_molecular_fraction(.2)==.1
 assert tumor_purity_to_molecular_fraction(1)==.5
 for x in [-.01,1.01,float('nan')]:
  with pytest.raises(ModelError):tumor_purity_to_molecular_fraction(x)

def test_calibration_manifest_rejects_test_samples(tmp_path):
 rs=[dict(donor=str(i),role='test',label='0',result='unused.json') for i in range(20)]
 with pytest.raises(ModelError,match='declared normal'):make_calibration(rs,tmp_path/'bank.json')

def test_reference_subset_removes_private_allele():
 panel={'donors':['a','a','b','b'],'hap_ids':['a|1','a|2','b|1','b|2'],'populations':['EAS']*4,
 'snps':np.zeros((2,4),dtype=int),'str_haps':{'L':np.array([0,2,0,1])},'locus_order':['L'],
 'loci':{'L':{'alleles':[{'repeat_sequence':'A'*n} for n in [10,11,12]]}},'phase_compatible':{'L':{'0':'ALL'}}}
 q=subset_panel(panel,['b'])
 assert q['donors']==['b','b'] and len(q['loci']['L']['alleles'])==2
 assert 'A'*12 not in [a['repeat_sequence'] for a in q['loci']['L']['alleles']]
