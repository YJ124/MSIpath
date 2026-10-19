"""Portable MSIpath command-line interface; the released model is the default."""
import argparse,json
from pathlib import Path
from . import __version__
from .model.common import ModelError,read_table
from .pipeline import infer,make_calibration

def main():
    parser=argparse.ArgumentParser(description='MSIpath: personalized reference and single-allele MSI evidence')
    parser.add_argument('--version',action='version',version='MSIpath '+__version__)
    sub=parser.add_subparsers(dest='command',required=True)
    call=sub.add_parser('call',help='BAM/CRAM/SAM or saved evidence to sample score and calibrated call')
    group=call.add_mutually_exclusive_group(required=True)
    group.add_argument('--bam',type=Path);group.add_argument('--evidence-dir',type=Path)
    for name in ['panel','normal-profile','output']:call.add_argument('--'+name,type=Path,required=True)
    call.add_argument('--reference',type=Path);call.add_argument('--loci',type=Path)
    call.add_argument('--calibration',type=Path)
    call.add_argument('--sample-id',required=True);call.add_argument('--donor-id',required=True)
    cal=sub.add_parser('calibrate',help='Build a bank from independently scored normal donors')
    cal.add_argument('--manifest',required=True,type=Path);cal.add_argument('--output',required=True,type=Path)
    ref=sub.add_parser('build-reference',help='Build schema-v2 panel from jointly phased SNP/STR data')
    ref.add_argument('--config',required=True,type=Path);ref.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    try:
        if args.command=='call':
            result=infer(panel_dir=args.panel,profile_file=args.normal_profile,output=args.output,
                sample_id=args.sample_id,donor_id=args.donor_id,bam=args.bam,reference=args.reference,
                evidence_dir=args.evidence_dir,loci_file=args.loci,calibration=args.calibration)
        elif args.command=='calibrate':
            records=read_table(args.manifest)
            for r in records:
                p=Path(r['result']);r['result']=str(p if p.is_absolute() else args.manifest.resolve().parent/p)
            result=make_calibration(records,args.output)
        else:
            import yaml
            from .model.reference import build_panel
            cfg=yaml.safe_load(args.config.read_text());base=args.config.resolve().parent
            for key in ['sample_metadata','loci','ancestry_markers','reference_fasta']:
                if cfg.get(key):
                    p=Path(cfg[key]);cfg[key]=str(p if p.is_absolute() else base/p)
            for key in ['gfa','registries','joint_phased_vcfs']:
                cfg[key]=[str(Path(p) if Path(p).is_absolute() else base/p) for p in cfg[key]]
            panel=build_panel(cfg,args.output);result=dict(panel=str(args.output),loci=len(panel['loci']))
        print(json.dumps(result,ensure_ascii=False,allow_nan=False));return 0
    except (ModelError,FileNotFoundError,KeyError,ValueError) as e:
        parser.exit(2,'MSIpath: '+str(e)+'\n')

if __name__=='__main__':raise SystemExit(main())
