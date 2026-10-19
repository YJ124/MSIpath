import argparse,json
from pathlib import Path
from msipath.pipeline import infer
ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
metadata=json.loads((ROOT/'demo_manifest.json').read_text())
for r in metadata:
 result=infer(panel_dir=ROOT/'research_panel',profile_file=ROOT/'normal_profiles.json',
     output=a.output/r['sample'],sample_id=r['sample'],donor_id=r['donor'],
     evidence_dir=ROOT/'evidence'/r['sample'],loci_file=ROOT/'common_loci.json',
     calibration=ROOT/'calibration_common15.json')
 assert abs(result['score']-r['expected_common15_score'])<1e-6
 print(r['sample'],result['call'],result['score'])
