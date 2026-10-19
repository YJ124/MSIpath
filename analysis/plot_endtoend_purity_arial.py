"""Redraw completed purity supplement; preserve all frozen analysis inputs."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
import numpy as np

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--font-dir', type=Path, help='Directory containing licensed Arial TTF files')
args = parser.parse_args()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'paper/figures/end_to_end_purity_grid'
OUT.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT / 'paper/data/end_to_end_purity_grid/performance_metrics.tsv'
SCORES = ROOT / 'paper/data/end_to_end_purity_grid/combined_sample_scores.tsv'
if args.font_dir:
    for name in ['arial.ttf', 'arialbd.ttf', 'ariali.ttf', 'arialbi.ttf']:
        font_manager.fontManager.addfont(str(args.font_dir / name))
font_manager.findfont(font_manager.FontProperties(family='Arial'), fallback_to_default=False)
plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 12, 'axes.labelsize': 13,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'axes.linewidth': 0.8, 'axes.unicode_minus': False,
})
with SOURCE.open() as f:
    metrics = list(csv.DictReader(f, delimiter='\t'))
with SCORES.open() as f:
    samples = list(csv.DictReader(f, delimiter='\t'))
METHODS = ['MSIpath', 'MSIsensor-pro', 'MSIsensor2', 'mSINGS', 'MANTIS2']
COLORS = dict(zip(METHODS, ['#D55E00', '#0072B2', '#E69F00', '#009E73', '#CC79A7']))
LABELS = {m: m for m in METHODS}
LABELS['MANTIS2'] = 'MANTIS2 (paired normal)'
PURITIES = [.2, .4, .6, .8]
PANELS = [('A', 'panel200', .1), ('B', 'panel200', .3),
          ('C', 'model_common', .1), ('D', 'model_common', .3)]
selected = []
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
for ax, (letter, scope, burden) in zip(axes.flat, PANELS):
    methods = [m for m in METHODS if scope == 'model_common' or m != 'MSIsensor2']
    width = .80 / len(methods)
    for j, method in enumerate(methods):
        rows = []
        for p in PURITIES:
            found = [r for r in metrics if r['cohort'] == 'all'
                     and r['scope'] == scope and r['method'] == method
                     and float(r['b']) == burden and float(r['p']) == p
                     and r['endpoint'] == 'sensitivity']
            assert len(found) == 1, (letter, method, p, len(found))
            row = found[0]
            # Independent agreement check against every individual sample call.
            raw = [r for r in samples if r['role'] == 'test'
                   and r['scope'] == scope and r['method'] == method
                   and float(r['b']) == burden and float(r['p']) == p]
            assert len(raw) == 200, (letter, method, p, len(raw))
            assert len({r['donor'] for r in raw}) == 200
            positives = sum(r['matched_call'] == 'MSI' for r in raw)
            assert int(row['instances']) == 200
            assert int(row['positives']) == positives
            assert abs(float(row['rate']) - positives / 200) < 1e-12
            rows.append(row)
            selected.append({'panel': letter, **row})
        y = np.array([float(r['rate']) * 100 for r in rows])
        lower = np.array([float(r['CI_low']) * 100 for r in rows])
        upper = np.array([float(r['CI_high']) * 100 for r in rows])
        x = np.arange(4) + (j - (len(methods) - 1) / 2) * width
        ax.bar(x, y, width=width * .91, color=COLORS[method], zorder=3,
               yerr=np.maximum(0, np.array([y - lower, upper - y])),
               error_kw={'ecolor': '#555555', 'elinewidth': .85,
                         'capsize': 2.1, 'capthick': .85})
    ax.set_xticks(np.arange(4), ['20', '40', '60', '80'])
    ax.set_xlabel('Tumor purity (%)')
    ax.set_ylabel('Sensitivity (%)')
    ax.set_ylim(0, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlim(-.57, 3.57)
    ax.grid(axis='y', color='#E8EBEE', linewidth=.8, zorder=0)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['bottom', 'left']].set_color('#555555')
    ax.tick_params(length=3, width=.7)
    ax.text(-.105, 1.035, letter, transform=ax.transAxes,
            fontsize=19, fontweight='bold', va='bottom')
fig.subplots_adjust(left=.075, right=.987, top=.95, bottom=.13,
                    hspace=.33, wspace=.22)
fig.legend(handles=[Patch(color=COLORS[m], label=LABELS[m]) for m in METHODS],
           loc='lower center', bbox_to_anchor=(.53, .033), ncol=5,
           frameon=False, fontsize=11.5, handlelength=1.2,
           handletextpad=.45, columnspacing=1.6)
base = OUT / 'end_to_end_purity_20_40_60_80_Arial'
fig.savefig(base.with_suffix('.png'), dpi=300, facecolor='white')
fig.savefig(base.with_suffix('.pdf'), facecolor='white',
            metadata={'Title': 'MSIpath end-to-end detection across four tumor purities',
                      'Subject': '200 donors per condition; Wilson 95% confidence intervals'})
plt.close(fig)
with (OUT / 'plotted_values.tsv').open('w') as f:
    writer = csv.DictWriter(f, fieldnames=list(selected[0]), delimiter='\t')
    writer.writeheader()
    writer.writerows(selected)
manifest = {
    'metric_source': str(SOURCE),
    'metric_source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'sample_source_sha256': hashlib.sha256(SCORES.read_bytes()).hexdigest(),
    'font_regular': font_manager.findfont(font_manager.FontProperties(family='Arial')),
    'font_bold': font_manager.findfont(font_manager.FontProperties(family='Arial', weight='bold')),
    'verified_bars': len(selected), 'denominator_per_bar': 200,
    'purity_grid': PURITIES,
    'panels': [{'panel': a, 'scope': b, 'burden': c} for a, b, c in PANELS],
    'error_bars': 'Wilson 95% confidence intervals',
    'unresolved_samples': 'Retained in the denominator; not counted as detected MSI',
}
(OUT / 'plot_manifest.json').write_text(json.dumps(manifest, indent=2))
(OUT / 'caption.txt').write_text(
    'End-to-end sample-level sensitivity at 30x sequencing depth. '
    'A: 200 STR loci, b=10%; B: 200 STR loci, b=30%; '
    'C: 15 common STR loci, b=10%; D: 15 common STR loci, b=30%. '
    'b is the proportion of altered STR loci. Each condition includes the same 200 test donors. '
    'Bars indicate the fraction classified as MSI using frozen independently calibrated thresholds. '
    'Unresolved samples remain in the denominator. Error bars are Wilson 95% confidence intervals. '
    'MSIsensor2 is evaluated in the common-locus panels C and D. '
    'Tumor purity follows the existing diploid, clonal, heterozygous simulation convention.\n')
print(json.dumps(manifest, indent=2))
print('Created:', base.with_suffix('.png'), base.with_suffix('.pdf'))
