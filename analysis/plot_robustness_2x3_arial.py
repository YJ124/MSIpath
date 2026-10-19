#!/usr/bin/env python3
"""Redraw Experiment 3: 2 rows x 3 columns, Arial, measured 60% purity.

Run: python plot_robustness_2x3_arial.py
Requires numpy, pandas, matplotlib, Pillow. No interpolation is used.
The 30x/60% point uses the supplement's exact 50 donor manifest, frozen
end-to-end scores, and frozen same-depth normal calibration (not the
200-donor aggregate). Other panels retain the original experiments.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image

PROJECT = Path(__file__).resolve().parents[1]
METHODS = ['MSIpath', 'MSIsensor-pro', 'mSINGS', 'MANTIS2']
COLORS = ['#D55E00', '#0072B2', '#009E73', '#CC79A7']
MARKERS = ['o', '^', 's', 'D']
STYLES = ['-', '--', ':', '-.']
DEPTHS = [10, 20, 30, 60]
PURITIES = [.1, .2, .4, .6, .8]


def read(path):
    return pd.read_csv(path, sep='\t')


def one(df, **filters):
    keep = pd.Series(True, index=df.index)
    for key, value in filters.items():
        keep &= np.isclose(df[key], value) if isinstance(value, float) else df[key].eq(value)
    out = df.loc[keep]
    assert len(out) == 1, (filters, len(out))
    return out.iloc[0]


def wilson(k, n):
    z = 1.95996398454
    p = k/n
    den = 1+z*z/n
    mid = (p+z*z/(2*n))/den
    half = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0., mid-half), min(1., mid+half)


def add_p60(data, supplement, output):
    paths = {
        'original_metrics': data/'robustness/performance_metrics.tsv',
        'original_scores': data/'robustness/sample_scores.tsv',
        'end_to_end_scores': data/'end_to_end/sample_scores.tsv',
        'supplement_manifest': supplement/'simulation/samples.tsv',
        'supplement_calibration': supplement/'config/calibration.json',
    }
    df = read(paths['original_metrics'])
    old = read(paths['original_scores'])
    source = read(paths['end_to_end_scores'])
    manifest = read(paths['supplement_manifest'])
    manifest = manifest.loc[manifest.depth.eq(30)]
    assert len(manifest) == 50 and manifest.donor.nunique() == 50
    assert np.allclose(manifest.p, .6) and np.allclose(manifest.b, .3)
    donors = set(manifest.donor)
    old_donors = set(old.loc[old.role.eq('test') & old.depth.eq(30) & old.group.eq('depth_fraction_grid'), 'donor'])
    assert donors == old_donors, 'Supplement must use the same 50 robustness donors'
    cal = pd.DataFrame(json.loads(paths['supplement_calibration'].read_text(encoding='utf-8')))
    records, sample_records = [], []
    for method in METHODS:
        ss = source.loc[source.method.eq(method) & source.scope.eq('panel200') &
                        source.role.eq('test') & source.depth.eq(30) &
                        np.isclose(source.b, .3) & np.isclose(source.p, .6) & source.donor.isin(donors)]
        assert len(ss) == 50 and set(ss.donor) == donors
        assert set(ss['sample']) == set(manifest['reuse'].str.replace('_F30', '_P60', regex=False))
        null = cal.loc[cal.method.eq(method) & cal.scope.eq('panel200') & cal.depth.eq(30) & cal.score.notna()]
        assert len(null) == 40 and not donors.intersection(null.donor)
        for r in ss.itertuples():
            assert pd.notna(r.score) and r.matched_calibration_n == 40
            p = (1+np.count_nonzero(null.score.to_numpy() >= r.score))/(len(null)+1)
            assert np.isclose(p, r.matched_p), (method, r.sample, p, r.matched_p)
            call = 'MSI' if p <= .05 else 'MSS'
            assert call == r.matched_call
            sample_records.append(dict(method=method, donor=r.donor, source_sample=r.sample,
                                       depth=30, tumor_purity=.6, altered_locus_fraction=.3,
                                       score=r.score, calibration_n=40, calibrated_p=p, call=call))
        k = int(ss.matched_call.eq('MSI').sum())
        lo, hi = wilson(k, 50)
        records.append(dict(method=method, scope='panel200', group='depth_fraction_grid',
                            condition='P60', depth=30, b=.3, p=.6, direction='mixed',
                            instances=50, endpoint='sensitivity', matched_rate=k/50,
                            matched_CI_low=lo, matched_CI_high=hi, matched_unresolved=0.))
    assert not ((df.scope == 'panel200') & (df.depth == 30) & np.isclose(df.p, .6) & (df.group == 'depth_fraction_grid')).any()
    extra = pd.DataFrame(records)
    pd.DataFrame(sample_records).to_csv(output/'Verified_purity60_sample_results.tsv', sep='\t', index=False)
    extra.to_csv(output/'Verified_purity60_metrics.tsv', sep='\t', index=False)
    (output/'Source_checksums.json').write_text(json.dumps({k: {'path': str(v), 'sha256': hashlib.sha256(v.read_bytes()).hexdigest()} for k, v in paths.items()}, indent=2), encoding='utf-8')
    print(extra[['method', 'instances', 'matched_rate', 'matched_CI_low', 'matched_CI_high']].to_string(index=False))
    return pd.concat([df, extra], ignore_index=True)


def values(rows):
    y, lo, hi = [], [], []
    for r in rows:
        if r.matched_unresolved == 1 or pd.isna(r.matched_rate):
            y.append(np.nan); lo.append(np.nan); hi.append(np.nan)
        else:
            v, l, h = r.matched_rate, r.matched_CI_low, r.matched_CI_high
            assert 0 <= l <= v+1e-12 and v <= h+1e-12 and h <= 1+1e-12
            y.append(v*100); lo.append(max(0, (v-l)*100)); hi.append(max(0, (h-v)*100))
    return y, np.array([lo, hi])


def curve(ax, df, queries, x, xlabel, ylabel='Sensitivity (%)', ymax=105):
    for method, color, marker, line in zip(METHODS, COLORS, MARKERS, STYLES):
        y, error = values([one(df, method=method, scope='panel200', **q) for q in queries])
        ax.errorbar(x, y, yerr=error, color=color, marker=marker, linestyle=line,
                    linewidth=1.7, markersize=4.8, elinewidth=.8, capsize=2, capthick=.8)
    ax.set(xlabel=xlabel, ylabel=ylabel, xticks=x, ylim=(0, ymax))
    ax.margins(x=.07)


def plot(df):
    for font in ['arial.ttf', 'arialbd.ttf', 'ariali.ttf', 'arialbi.ttf']:
        path = Path('C:/Windows/Fonts')/font
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    font = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 11, 'axes.labelsize': 11,
                         'xtick.labelsize': 10, 'ytick.labelsize': 10, 'axes.linewidth': .8,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'legend.fontsize': 11, 'legend.frameon': False,
                         'figure.facecolor': 'white', 'savefig.facecolor': 'white',
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 3, figsize=(17.4, 8.2), layout='constrained')
    a, b, c, d, e, f = axes.flat
    grid = 'depth_fraction_grid'
    curve(a, df, [dict(group=grid, depth=x, p=.2) for x in DEPTHS], DEPTHS, 'Sequencing depth (fold coverage)')
    curve(b, df, [dict(group=grid, depth=30, p=x) for x in PURITIES], [int(x*100) for x in PURITIES], 'Tumor purity (%)')
    curve(c, df, [dict(group=grid, depth=x, condition='N') for x in DEPTHS], DEPTHS,
          'Sequencing depth (fold coverage)', 'False-positive rate (%)', 35)
    for method, color, marker, line in zip(METHODS, COLORS, MARKERS, STYLES):
        y = [100*one(df, method=method, scope='panel200', group=grid, depth=x, condition='N').matched_unresolved for x in DEPTHS]
        d.plot(DEPTHS, y, color=color, marker=marker, linestyle=line, linewidth=1.7, markersize=4.8)
    d.set(xlabel='Sequencing depth (fold coverage)', ylabel='Unresolved normal samples (%)',
          xticks=DEPTHS, yticks=[0, 25, 50, 75, 100], ylim=(-2, 105))
    d.margins(x=.07)
    for ax, p in [(e, .2), (f, .8)]:
        queries = [dict(group=grid, depth=30, p=p),
                   dict(group='direction_control', depth=30, p=p, condition='expansion'),
                   dict(group='direction_control', depth=30, p=p, condition='contraction')]
        for j, (method, color) in enumerate(zip(METHODS, COLORS)):
            y, error = values([one(df, method=method, scope='panel200', **q) for q in queries])
            ax.bar(np.arange(3)+(j-1.5)*.195, y, width=.1794, yerr=error,
                   color=color, edgecolor='white', linewidth=.5,
                   error_kw={'elinewidth': .9, 'capsize': 2, 'capthick': .9, 'ecolor': '#444444'})
        ax.set(xticks=range(3), xticklabels=['Mixed\ndirections', 'Expansion\nonly', 'Contraction\nonly'],
               xlabel='Direction of repeat-length alteration', ylabel='Sensitivity (%)',
               ylim=(0, 105), yticks=[0, 25, 50, 75, 100])
    for label, ax in zip('ABCDEF', axes.flat):
        ax.text(-.11, 1.025, label, transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom')
        ax.set_axisbelow(True)
        ax.grid(axis='y', color='#E4E7E8', linewidth=.65)
        ax.tick_params(length=3)
    fig.legend(handles=[Line2D([], [], color=color, marker=marker, linestyle=line,
                             linewidth=1.7, label='MANTIS2 (paired normal)' if method == 'MANTIS2' else method)
                        for method, color, marker, line in zip(METHODS, COLORS, MARKERS, STYLES)],
               loc='outside lower center', ncol=4, columnspacing=2.)
    fig.get_layout_engine().set(w_pad=.09, h_pad=.09, wspace=.035, hspace=.06)
    print('Verified Arial font:', font)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=PROJECT/'paper/data')
    parser.add_argument('--supplement-dir', type=Path, default=PROJECT/'paper/data/robustness_purity60')
    parser.add_argument('--output-dir', type=Path, default=PROJECT/'paper/output/figures/Robustness_2x3_Arial')
    parser.add_argument('--dpi', type=int, default=600)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = add_p60(args.data_dir, args.supplement_dir, args.output_dir)
    fig = plot(df)
    path = args.output_dir/'Experiment3_2x3_Arial_purity60.png'
    fig.savefig(path, dpi=args.dpi, bbox_inches='tight', pad_inches=.12)
    plt.close(fig)
    with Image.open(path) as im:
        im.thumbnail((1900, 1000))
        im.save(args.output_dir/'Preview.png')
    notes = '''# Panel definitions and provenance

A: Sensitivity versus depth; tumor purity 20%.
B: Sensitivity versus tumor purity (10%, 20%, 40%, 60%, 80%); depth 30x.
C: Ordinary-normal false-positive rate versus depth.
D: Unresolved ordinary-normal sample fraction versus depth.
E: Direction controls at 20% tumor purity, 30x depth.
F: Direction controls at 80% tumor purity, 30x depth.

All panels use 200 STR loci, 50 donors per condition; altered samples have
30% altered loci. Bars/points use all 50 samples as denominator; error bars
are 95% Wilson intervals. Fully unresolved conditions are gaps, not zeros.
Each method uses its frozen same-depth independent normal calibration.
Panel letters are assigned left to right, top to bottom. No panel titles.

The added 60% purity point is derived from exactly the 50 donors specified
in the purity_p60_v1 supplement manifest. Their original B30_P60 scores
are reused as explicitly specified by that experiment. All 200 donor-method
records are independently recalibrated against the supplement's 40 frozen
normal scores and agree with released per-sample p-values and calls.
No 200-donor aggregate and no interpolation is used.
The other-depth supplement is incomplete in the local status snapshot and
is not used. Pure expansion/contraction at 60% was not added to E or F.
'''
    (args.output_dir/'Panel_definitions.md').write_text(notes, encoding='utf-8')
    print(path)


if __name__ == '__main__':
    main()
