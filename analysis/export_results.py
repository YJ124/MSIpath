"""Import frozen experiments, rename the released model, and convert simulation purity."""
import argparse,csv,hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
METHODS={'MSIpath-B':'MSIpath','MSIsensor-pro':'MSIsensor-pro','MSIsensor2':'MSIsensor2','mSINGS':'mSINGS','MANTIS2_paired':'MANTIS2'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))
def label(s):return re.sub(r'(^|_)F(\d+)(?=_|$)',lambda m:m[1]+'P'+f'{int(m[2])*2:02}',s)
def transform(r):
 out={}
 for k,v in r.items():
  if k in ['evidence_dir','bh_significant']:continue
  if k=='method':out[k]=METHODS[v]
  elif k=='f':out['p']=format(2*float(v),'.12g')
  elif k=='burden':out['b']=v
  elif k in ['sample','paired_normal','condition']:out[k]=label(v)
  else:out[k]=v
 return out
def table(p,data):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',encoding='utf8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
def main():
 p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);a=p.parse_args();base=a.workspace; audit=[]
 for source,target in [('experiment','end_to_end'),('experiment3','mechanism'),('experiment4','robustness')]:
  for name in ['performance_metrics.tsv','sample_scores.tsv']:
   f=base/source/'B_v1/reports'/name;rs=read(f);rs=[transform(r) for r in rs if r['method'] in METHODS]
   assert all('f' not in r and r['method'] in METHODS.values() for r in rs)
   table(ROOT/'paper/data'/target/name,rs)
   audit.append(dict(source=source+'/reports/'+name,sha256=sha(f),exported_rows=len(rs)))
 res=base/'experiment2/B_v1/unified_resources'
 assert json.loads((res/'reports/COMPLETE.json').read_text())['executions']==270
 for name in ['resource_comparison.tsv','measurements.tsv']:
  f=res/'reports'/name;rs=[transform(r) for r in read(f)];table(ROOT/'paper/data/resources'/name,rs)
  audit.append(dict(source='unified_resources/reports/'+name,sha256=sha(f),exported_rows=len(rs)))
 hardware=json.loads((res/'reports/hardware.json').read_text())
 hardware['uname'][1]='local-benchmark-host'
 (ROOT/'paper/data/resources/hardware.json').write_text(json.dumps(hardware,indent=2)+'\n')
 stages=read(base/'experiment2/B_v1/reports/resource_summary.tsv')
 for r in stages:r['stage']=r['stage'].replace('Step4_B_','Step4_').replace('Step5_B_','Step5_')
 table(ROOT/'paper/data/resources/stage_scaling.tsv',stages)
 # Preserve within-MSIpath mechanism diagnostics in a separate appendix table.
 f=base/'experiment3/B_v1/reports/mechanism_locus_metrics.tsv';rs=read(f)
 modes={'MSIpath-B':'personalized_mixture','B_MAP':'single_MAP_genotype','B_population_full_panel':'full_panel_population'}
 for r in rs:r['mode']=modes[r['mode']]
 table(ROOT/'paper/data/mechanism/within_method_diagnostics.tsv',rs)
 (ROOT/'paper/data/provenance.json').write_text(json.dumps(dict(files=audit,method_mapping=METHODS,
    parameter_conversion='tumor purity = 2 * source abnormal molecule fraction; clonal heterozygous, diploid copy-neutral simulation',
    no_scores_or_calls_changed=True,original_results_preserved=True),indent=2)+'\n')
 print('Exported final five-method data and unified resources')
if __name__=='__main__':main()
