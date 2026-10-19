#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MSIpath publication plots and Markdown tables (revision 2).

Dependencies: Python >= 3.10, numpy, pandas, matplotlib >= 3.8.
Run:
    python plot_three_experiments.py
    python plot_three_experiments.py --dpi 600 --formats png pdf
    python plot_three_experiments.py --no-heatmap

Outputs: three main figures; an optional supplementary heatmap; three Markdown
results tables and panel descriptions. No resource plots, resource tables or
PowerPoint files are produced. There are no figure titles, point-value labels
or significance symbols. Panel letters, full axis names and legends remain.

Data conventions:
* Released tumor purity equals twice altered-molecule fraction only under the
  diploid, copy-neutral, clonal heterozygous simulation assumption.
* Sensitivity and false-positive rates use same-depth independent-normal
  calibration, not identical realized specificity. Wilson intervals include
  all test instances in the denominator. Fully unresolved results remain NaN
  in the plots; they are never substituted with valid zero detection.
* Ranking metrics use scored samples only. The AUPRC field stores average
  precision, not trapezoidal integration of an interpolated precision-recall
  curve. Bootstrap intervals are copied from the frozen input tables.
* Rare and novel inherited-normal edits retain the flanking SNP background:
  these engineered stress tests do not estimate population-wide false positives.
* Full-panel normal modelling changes allele support and weights jointly.
* Exact STR fragments count filtered paired-fragment observations, not unique
  original molecules. Their sample identities and callable-locus counts are
  cross-checked against the released per-sample scores before plotting.
"""
from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
METHODS = ["MSIpath", "MSIsensor-pro", "mSINGS", "MANTIS2"]
ALL_METHODS = ["MSIpath", "MSIsensor-pro", "MSIsensor2", "mSINGS", "MANTIS2"]
DISPLAY = {"MANTIS2": "MANTIS2 (paired normal)"}
COLORS = {"MSIpath": "#D55E00", "MSIsensor-pro": "#0072B2",
          "MSIsensor2": "#E69F00", "mSINGS": "#009E73", "MANTIS2": "#CC79A7"}
MARKERS = dict(zip(ALL_METHODS, ["o", "^", "P", "s", "D"]))
LINE_STYLES = {"MSIpath": "-", "MSIsensor-pro": "--", "MSIsensor2": "-",
               "mSINGS": ":", "MANTIS2": "-."}
GROUPS = ["N", "allele_imbalance", "rare_in_panel", "novel_normal", "novel_somatic"]
GROUP_NAMES = [
    "Ordinary normal background",
    "Allelic imbalance (80:20)",
    "Rare inherited alleles within the reference panel",
    "Novel inherited alleles outside the reference panel",
    "Novel somatic alterations",
]
GROUP_TICKS = [
    "Ordinary normal\nbackground", "Allelic imbalance\n(80:20)",
    "Rare inherited alleles\nwithin the reference panel",
    "Novel inherited alleles\noutside the reference panel",
    "Novel somatic\nalterations",
]
GROUP_COLORS = ["#6C757D", "#0072B2", "#E69F00", "#9467BD", "#D55E00"]
MODES = {
    "personalized_mixture": "Personalized normal genotype mixture",
    "single_MAP_genotype": "Maximum a posteriori normal genotype",
    "full_panel_population": "Full reference panel population model",
}
DEPTHS = [10, 20, 30, 60]
PURITIES = [.1, .2, .4, .8]
CONDITIONS = ["B10_P20", "B10_P60", "B30_P20", "B30_P60"]


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelsize": 11, "axes.linewidth": .8,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10, "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "pdf.fonttype": 42,
        "ps.fonttype": 42, "svg.fonttype": "none",
    })


def read(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def select(df, **filters):
    mask = pd.Series(True, index=df.index)
    for key, value in filters.items():
        mask &= (np.isclose(df[key], value) if isinstance(value, float) else df[key].eq(value))
    return df.loc[mask]


def one(df, **filters):
    result = select(df, **filters)
    if len(result) != 1:
        raise ValueError(f"Expected exactly one record: {filters}; found {len(result)}")
    return result.iloc[0]


def axes_style(ax, letter, grid=True):
    # Letters are the only panel annotations; no titles or numerical annotations.
    ax.text(-.10, 1.025, letter, transform=ax.transAxes, fontsize=16,
            fontweight="bold", ha="left", va="bottom", clip_on=False)
    ax.set_axisbelow(True)
    if grid:
        ax.grid(axis="y", color="#E4E7E8", linewidth=.65)
    ax.tick_params(length=3)


def legend(fig, methods, bar=False):
    handles = ([Patch(facecolor=COLORS[m], label=DISPLAY.get(m, m)) for m in methods]
               if bar else [Line2D([], [], color=COLORS[m], marker=MARKERS[m], linestyle=LINE_STYLES[m],
                                   linewidth=1.8, label=DISPLAY.get(m, m)) for m in methods])
    fig.legend(handles=handles,
               loc="outside lower center", ncol=len(methods), columnspacing=1.7)


def rate_values(rows):
    estimates, lower, upper = [], [], []
    for r in rows:
        if r["matched_unresolved"] == 1 or pd.isna(r["matched_rate"]):
            estimates.append(np.nan); lower.append(np.nan); upper.append(np.nan)
            continue
        v, lo, hi = [float(r[k]) for k in ("matched_rate", "matched_CI_low", "matched_CI_high")]
        if not -1e-10 <= lo <= v+1e-10 <= hi+1e-10 <= 1+1e-10:
            raise ValueError("Invalid rate confidence interval")
        estimates.append(v*100); lower.append(max(0, (v-lo)*100)); upper.append(max(0, (hi-v)*100))
    return np.array(estimates), np.array([lower, upper])


def bars(ax, df, queries, methods, labels, xlabel, ylabel="Sensitivity (%)"):
    width = .78/len(methods)
    for index, method in enumerate(methods):
        values, errors = rate_values([one(df, method=method, **q) for q in queries])
        ax.bar(np.arange(len(labels)) + (index-(len(methods)-1)/2)*width,
               values, width=width*.92, yerr=errors,
               color=COLORS[method], edgecolor="white", linewidth=.5,
               error_kw={"elinewidth": .9, "capsize": 2, "capthick": .9,
                         "ecolor": "#444444"})
    ax.set(xticks=np.arange(len(labels)), xticklabels=labels,
           xlabel=xlabel, ylabel=ylabel, ylim=(0, 105), yticks=[0, 25, 50, 75, 100])


def curves(ax, df, queries, methods, x, xlabel, ylabel="Sensitivity (%)", ymax=105):
    for method in methods:
        values, errors = rate_values([one(df, method=method, **q) for q in queries])
        # NaNs break lines. Missing/fully unresolved points do not become zeros.
        ax.errorbar(x, values, yerr=errors, color=COLORS[method], marker=MARKERS[method], linestyle=LINE_STYLES[method],
                    linewidth=1.8, markersize=5, elinewidth=.8, capsize=2, capthick=.8)
    ax.set(xlabel=xlabel, ylabel=ylabel, xticks=x, ylim=(0, ymax))
    ax.margins(x=.07)


def boxes(ax, arrays, labels, colors, seed=20261013):
    rng = np.random.default_rng(seed)
    for index, (values, color) in enumerate(zip(arrays, colors)):
        values = np.asarray(values, dtype=float)
        if not len(values) or not np.isfinite(values).all():
            raise ValueError("Box plots require complete observed scores")
        ax.boxplot([values], positions=[index], widths=.48, whis=1.5,
                   showfliers=False, patch_artist=True,
                   boxprops={"facecolor": "white", "edgecolor": color, "linewidth": 1.5},
                   medianprops={"color": color, "linewidth": 2},
                   whiskerprops={"color": color, "linewidth": 1.1},
                   capprops={"color": color, "linewidth": 1.1})
        # All observations, including outliers, remain visible.
        ax.scatter(index+rng.uniform(-.15, .15, len(values)), values, color=color,
                   s=13, alpha=.65, linewidths=0, zorder=3)
    ax.set(xticks=range(len(labels)), xticklabels=labels, xlim=(-.6, len(labels)-.4))


def experiment1(df):
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), layout="constrained")
    for index, (scope, methods, burden) in enumerate([
        ("panel200", METHODS, 10), ("panel200", METHODS, 30),
        ("model_common", ALL_METHODS, 10), ("model_common", ALL_METHODS, 30),
    ]):
        ax = axes.flat[index]
        queries = [{"scope": scope, "condition": f"B{burden}_P{p}"} for p in [20, 60]]
        bars(ax, df, queries, methods, ["20", "60"], "Tumor purity (%)")
        axes_style(ax, chr(65+index))
    legend(fig, ALL_METHODS, bar=True)
    return fig


def experiment2(scores):
    fig, axes = plt.subplots(2, 2, figsize=(15.4, 10.8), layout="constrained")
    for index, (ax, method) in enumerate(zip(axes.flat, METHODS)):
        arrays = []
        for condition in GROUPS:
            rr = select(scores, role="test", scope="panel200", method=method, condition=condition)
            if len(rr) != 50 or rr["sample"].duplicated().any():
                raise ValueError(f"Mechanism group must contain 50 distinct samples: {method}, {condition}")
            arrays.append(rr.sort_values("sample")["score"].to_numpy())
        boxes(ax, arrays, GROUP_TICKS, GROUP_COLORS)
        ax.set(xlabel="Simulated biological condition", ylabel="Sample-level evidence score")
        # Method identity is a legend, not a figure or panel title.
        ax.legend(handles=[Patch(facecolor="none", edgecolor="none", label=DISPLAY.get(method, method))],
                  loc="upper left", handlelength=0, handletextpad=0)
        plt.setp(ax.get_xticklabels(), rotation=20, ha="right", rotation_mode="anchor")
        axes_style(ax, chr(65+index))
        ax.margins(y=.14)
    return fig


def release_id(value):
    return re.sub(r"(^|_)F(\d+)(?=_|$)", lambda m: m[1]+"P"+f"{int(m[2])*2:02d}", str(value))


def fragment_data(path, released):
    raw = read(path).copy()
    for col in ["sample", "condition"]:
        raw[col] = raw[col].map(release_id)
    raw["p"] = 2*raw["f"] if "f" in raw else raw["p"]
    chosen = select(raw, condition="P20", p=.2)
    published = select(released, role="test", method="MSIpath", scope="panel200",
                       group="depth_fraction_grid", p=.2)
    joined = chosen.merge(published[["sample", "donor", "depth", "callable_loci"]],
                          on=["sample", "donor", "depth"], suffixes=("_raw", "_released"),
                          how="outer", indicator=True, validate="one_to_one")
    if len(joined) != 200 or not joined["_merge"].eq("both").all():
        raise ValueError("Fragment data do not match released sample identities")
    if not np.array_equal(joined["callable_loci_raw"], joined["callable_loci_released"]):
        raise ValueError("Fragment data do not match released callable-locus counts")
    return chosen


def direction_queries(p):
    return [
        {"scope": "panel200", "group": "depth_fraction_grid", "depth": 30, "p": p},
        {"scope": "panel200", "group": "direction_control", "depth": 30, "p": p, "condition": "expansion"},
        {"scope": "panel200", "group": "direction_control", "depth": 30, "p": p, "condition": "contraction"},
    ]


def experiment3(df):
    fig, axes = plt.subplots(3, 2, figsize=(13.8, 13), layout="constrained")
    curves(axes[0, 0], df, [
        {"scope": "panel200", "group": "depth_fraction_grid", "depth": d, "p": .2} for d in DEPTHS
    ], METHODS, DEPTHS, "Sequencing depth (fold coverage)")
    curves(axes[0, 1], df, [
        {"scope": "panel200", "group": "depth_fraction_grid", "depth": 30, "p": p} for p in PURITIES
    ], METHODS, [p*100 for p in PURITIES], "Tumor purity (%)")
    curves(axes[1, 0], df, [
        {"scope": "panel200", "group": "depth_fraction_grid", "depth": d, "condition": "N"} for d in DEPTHS
    ], METHODS, DEPTHS, "Sequencing depth (fold coverage)", "False-positive rate (%)", 35)
    for method in METHODS:
        values = [100*one(df, method=method, scope="panel200", group="depth_fraction_grid",
                         depth=d, condition="N")["matched_unresolved"] for d in DEPTHS]
        axes[1, 1].plot(DEPTHS, values, color=COLORS[method], marker=MARKERS[method], linestyle=LINE_STYLES[method],
                        linewidth=1.8, markersize=5)
    axes[1, 1].set(xlabel="Sequencing depth (fold coverage)", ylabel="Unresolved normal samples (%)",
                   xticks=DEPTHS, yticks=[0, 25, 50, 75, 100], ylim=(-2, 105))
    axes[1, 1].margins(x=.07)
    for ax, p in zip(axes[2], [.2, .8]):
        bars(ax, df, direction_queries(p), METHODS,
             ["Mixed directions", "Expansion only", "Contraction only"],
             "Direction of repeat-length alteration")
    for index, ax in enumerate(axes.flat):
        axes_style(ax, chr(65+index))
    legend(fig, METHODS)
    return fig


def supplementary_heatmap(df):
    # Plasma is perceptually uniform, with a strong purple-to-yellow lightness range.
    # Identical values have identical colours; no nonlinear exaggeration is used.
    cmap = plt.get_cmap("plasma").copy()
    cmap.set_bad("#B7BDC5")
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 9.2), layout="constrained")
    for index, (ax, method) in enumerate(zip(axes.flat, METHODS)):
        matrix = np.full((4, 4), np.nan)
        for row, p in enumerate(PURITIES[::-1]):
            for col, d in enumerate(DEPTHS):
                r = one(df, scope="panel200", method=method, group="depth_fraction_grid", depth=d, p=p)
                if r["matched_unresolved"] < 1:
                    matrix[row, col] = 100*r["matched_rate"]
        im = ax.imshow(np.ma.masked_invalid(matrix), cmap=cmap, vmin=0, vmax=100,
                       aspect="auto", interpolation="nearest")
        ax.set(xticks=range(4), xticklabels=DEPTHS,
               yticks=range(4), yticklabels=[int(p*100) for p in PURITIES[::-1]],
               xlabel="Sequencing depth (fold coverage)", ylabel="Tumor purity (%)")
        ax.set_xticks(np.arange(-.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-.5, 4, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.2)
        ax.tick_params(which="both", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
        axes_style(ax, chr(65+index), grid=False)
    fig.colorbar(im, ax=list(axes.flat), fraction=.035, pad=.025, label="Sensitivity (%)")
    fig.legend(handles=[Patch(facecolor="#B7BDC5", label="All samples unresolved")],
               loc="outside lower center")
    return fig


def fmt(x, digits=3):
    return "Not estimable" if pd.isna(x) else f"{float(x):.{digits}f}"


def interval(row, metric="matched_rate", lower="matched_CI_low", upper="matched_CI_high", scale=100, digits=1):
    if pd.isna(row[metric]):
        return "Not estimable"
    return f"{row[metric]*scale:.{digits}f} [{row[lower]*scale:.{digits}f}, {row[upper]*scale:.{digits}f}]"


def md_table(headers, rows):
    def escape(value):
        return str(value).replace("|", r"\|").replace("\n", " ")
    return "\n".join(["| "+" | ".join(map(escape, headers))+" |",
                      "| "+" | ".join(["---"]*len(headers))+" |"]+
                     ["| "+" | ".join(map(escape, row))+" |" for row in rows])+"\n\n"


def tables(out, primary, mechanism, diagnostic, robust, fragments):
    # Full metric names are deliberately retained in titles and column headings.
    text = "# Table 1. End-to-end detection performance and normal-sample classification\n\n"
    text += ("Rates use independent normal calibration at the matching sequencing depth (30-fold coverage). "
             "All rate denominators include unresolved samples. Intervals are 95% Wilson intervals. "
             "Ranking metrics use scored samples; average precision is the source AUPRC field. "
             "Correct normal classification is not defined as one minus the false-positive rate when unresolved calls occur.\n\n")
    for scope, methods, scope_name in [
        ("panel200", METHODS, "200 short tandem repeat loci"),
        ("model_common", ALL_METHODS, "15 commonly supported short tandem repeat loci"),
    ]:
        text += f"## {scope_name}: positive samples\n\n"
        rows = []
        for condition in CONDITIONS:
            for method in methods:
                r = one(primary, method=method, scope=scope, condition=condition)
                rows.append([DISPLAY.get(method, method), int(r["b"]*100), int(r["p"]*100),
                             int(r["instances"]), fmt(r["AUROC"]), fmt(r["AUPRC"]),
                             interval(r), fmt(r["matched_unresolved"]*100, 1),
                             fmt(r["score_available_rate"]*100, 1)])
        text += md_table(["Method", "Proportion of altered short tandem repeat loci (%)", "Tumor purity (%)", "Test instances",
                          "Area under the receiver operating characteristic curve",
                          "Average precision", "Sensitivity (%; 95% confidence interval)",
                          "Unresolved samples (%)", "Samples with ranking scores (%)"], rows)
        text += f"## {scope_name}: ordinary normal samples\n\n"
        rows = []
        for method in methods:
            r = one(primary, method=method, scope=scope, condition="N")
            rows.append([DISPLAY.get(method, method), int(r["instances"]), interval(r),
                         fmt(r["matched_correct_MSS_rate"]*100, 1), fmt(r["matched_unresolved"]*100, 1)])
        text += md_table(["Method", "Normal test instances",
                          "False-positive rate (%; 95% confidence interval)",
                          "Correctly classified normal samples (%)", "Unresolved samples (%)"], rows)
    (out/"Table1_End_to_end_performance.md").write_text(text, encoding="utf-8")

    text = "# Table 2. Normal-background stress tests and normal-model comparisons\n\n"
    text += ("Engineered rare and novel inherited-normal alleles retain the flanking single-nucleotide variant background. "
             "These challenges do not estimate false positives in an unselected population. "
             "Full-panel modelling changes supported alleles and their weights jointly.\n\n")
    text += "## Sample-level classification in 200 short tandem repeat loci\n\n"
    rows = []
    for group, name in zip(GROUPS, GROUP_NAMES):
        for method in METHODS:
            r = one(mechanism, method=method, scope="panel200", condition=group)
            rows.append([name, DISPLAY.get(method, method), int(r["instances"]),
                         interval(r), fmt(r["matched_unresolved"]*100, 1)])
    text += md_table(["Biological condition", "Method", "Test instances",
                      "Microsatellite-instability calls (%; 95% Wilson confidence interval)",
                      "Unresolved samples (%)"], rows)
    text += ("Calls in all inherited-normal conditions are false positives; calls in the somatic condition measure sensitivity.\n\n")
    for stratum, name in [("intervention", "Designated intervention loci"), ("background", "Unmodified background loci")]:
        text += f"## {name}: locus-level model comparison\n\n"
        rows = []
        for group, group_name in zip(["common_normal", *GROUPS[1:]], GROUP_NAMES):
            for mode, mode_name in MODES.items():
                r = one(diagnostic, mode=mode, group=group, stratum=stratum)
                rows.append([group_name, mode_name, int(r["donors"]), int(r["loci"]),
                             interval(r, "nominal_p05", "CI_low", "CI_high"),
                             fmt(r["uncalibrated_rate"]*100, 1)])
        text += md_table(["Biological condition", "Normal-background model", "Donors",
                          "Evaluated locus observations",
                          "Proportion of loci with nominal p-value ≤ 0.05 (%; 95% confidence interval)",
                          "Loci without calibration (%)"], rows)
    text += ("Locus fractions are averaged per donor. Intervals are the source 1,000-draw donor bootstrap "
             "95% intervals. Nominal locus significance is not the sample-level calling rule.\n")
    (out/"Table2_Mechanism_comparisons.md").write_text(text, encoding="utf-8")

    text = "# Table 3. Repeat-length direction comparisons and low-information detection\n\n"
    text += ("Direction comparisons use 200 short tandem repeat loci, 30-fold coverage, a 30% altered-locus proportion "
             "and 50 test instances per condition. MANTIS2 uses a paired normal sample. Receiver operating characteristic "
             "area intervals are the source donor-bootstrap intervals; sensitivity intervals are 95% Wilson intervals. "
             "Expansion and contraction share edited targets, haplotypes and absolute shift lengths; "
             "the mixed-direction condition does not always have the identical edited-locus set.\n\n")
    for purity in [.2, .8]:
        text += f"## Tumor purity of {int(purity*100)}%: direction-dependent discrimination\n\n"
        rows = []
        for query, direction in zip(direction_queries(purity), ["Mixed directions", "Expansion only", "Contraction only"]):
            for method in METHODS:
                r = one(robust, method=method, **query)
                rows.append([direction, DISPLAY.get(method, method),
                             interval(r, "AUROC", "AUROC_CI_low", "AUROC_CI_high", 1, 3),
                             fmt(r["AUPRC"]), interval(r), fmt(r["matched_unresolved"]*100, 1)])
        text += md_table(["Repeat-length alteration direction", "Method",
                          "Area under the receiver operating characteristic curve (95% confidence interval)",
                          "Average precision", "Sensitivity (%; 95% confidence interval)",
                          "Unresolved samples (%)"], rows)
    for method in METHODS:
        text += f"## {DISPLAY.get(method, method)}: sequencing depth and tumor purity\n\n"
        rows = []
        for p in PURITIES:
            cells = []
            for depth in DEPTHS:
                r = one(robust, method=method, scope="panel200", group="depth_fraction_grid", depth=depth, p=p)
                cells.append("All samples unresolved" if r["matched_unresolved"] == 1 else interval(r))
            rows.append([int(p*100), *cells])
        text += md_table(["Tumor purity (%)"]+[f"{d}-fold coverage: sensitivity (%; 95% confidence interval)" for d in DEPTHS], rows)
    text += ("An all-unresolved condition is not displayed as an observed zero sensitivity. "
             "Connecting segments in the figures are visual guides between measured conditions, not additional experiments.\n")
    text += "\n## Exact short tandem repeat fragment observations at 20% tumor purity\n\n"
    rows = []
    for depth in DEPTHS:
        values = select(fragments, depth=depth)["exact_STR_fragments"].to_numpy(dtype=float)
        if len(values) != 50 or not np.isfinite(values).all():
            raise ValueError("Expected 50 complete fragment observations per sequencing depth")
        rows.append([depth, len(values), fmt(np.median(values), 1),
                     fmt(np.quantile(values, .25), 1), fmt(np.quantile(values, .75), 1),
                     int(values.min()), int(values.max())])
    text += md_table(["Sequencing depth (fold coverage)", "Test instances",
                      "Median exact short tandem repeat fragments per sample",
                      "First quartile", "Third quartile", "Minimum", "Maximum"], rows)
    text += ("Quantiles use linear interpolation. Counts are filtered paired-fragment observations, "
             "not independent original molecules. Each sample identity and callable-locus count "
             "has been cross-checked against the released sample scores.\n")
    (out/"Table3_Direction_and_low_information.md").write_text(text, encoding="utf-8")


PANEL_GUIDE = """# Figure panel definitions

## Experiment 1

- A: 200 loci, 10% altered loci.
- B: 200 loci, 30% altered loci.
- C: 15 commonly supported loci, 10% altered loci.
- D: 15 commonly supported loci, 30% altered loci.

All panels compare 20% and 60% tumor purity at 30-fold coverage, with 200 test instances per condition.
Bars show sensitivity; error bars show 95% Wilson confidence intervals. Five methods appear in the common-locus analysis,
and four in the 200-locus analysis. Ranking metrics and normal outcomes are in Table 1.

## Experiment 2

- A: MSIpath.
- B: MSIsensor-pro.
- C: mSINGS.
- D: MANTIS2 with paired normal input.

Each group has 50 scored test samples at 30-fold coverage in the 200-locus scope. Axes retain method-specific raw score
scales. Boxes show the 25th–75th percentiles, central lines the median, and whiskers the most extreme observations within
1.5 interquartile ranges. Every observation is shown, including outliers; horizontal jitter is deterministic.
The 80:20 group changes the relative contributions of two alleles without introducing a somatic repeat-length change.
Rare inherited alleles are present in the full reference panel; novel inherited alleles are outside that panel.
These engineered conditions preserve the single-nucleotide variant background.

## Experiment 3

- A: Sensitivity versus sequencing depth at 20% tumor purity.
- B: Sensitivity versus tumor purity at 30-fold coverage.
- C: Ordinary-normal false-positive rate versus sequencing depth.
- D: Unresolved classification fraction in ordinary normal samples versus sequencing depth.
- E: Repeat-length direction comparison at 20% tumor purity and 30-fold coverage.
- F: Repeat-length direction comparison at 80% tumor purity and 30-fold coverage.

Positive conditions have 200 loci, 30% altered loci and 50 test instances. Sensitivity and false-positive rate uncertainty uses 95% Wilson intervals. Panel D shows observed fractions without confidence intervals.
Fully unresolved points are omitted; the corresponding table states their status explicitly. Fragments are filtered
paired-fragment observations, not independent original molecules. Their sample identities and callable-locus counts
match the released sample scores. No unseen experimental conditions are interpolated as measurements.

## Supplementary heatmap

A–D correspond to MSIpath, MSIsensor-pro, mSINGS and paired-normal MANTIS2, respectively.
All panels share a linear 0–100% sensitivity scale using the perceptually uniform Plasma colour map.
Grey means all samples unresolved and differs from the dark-purple colour for valid zero sensitivity.
Panels contain no numerical cell annotations.

## Scope and reproducibility

Tumor purity is twice the altered-molecule fraction only under the simulated diploid, copy-neutral, clonal heterozygous
alteration assumption. Error bars do not establish independence between perturbations of the same donor.
Rare/novel inherited-normal challenges remain limitations; the figures do not assert universal method superiority.
The source contains no resource-comparison plots or tables and creates no PowerPoint files.
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path, default=PROJECT/"paper/data")
    parser.add_argument("--fragment-data", type=Path,
                        default=PROJECT/"paper/data/robustness/effective_STR_fragments.tsv")
    parser.add_argument("--out-dir", type=Path, default=PROJECT/"paper/output/figures/Python_publication_revision")
    parser.add_argument("--formats", nargs="+", choices=["png", "pdf", "svg"], default=["png"])
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--no-heatmap", action="store_true")
    args = parser.parse_args()
    if args.dpi < 1:
        parser.error("DPI must be positive")
    style()
    source_files = {
        "primary": args.data_dir/"end_to_end/performance_metrics.tsv",
        "mechanism": args.data_dir/"mechanism/performance_metrics.tsv",
        "mechanism_scores": args.data_dir/"mechanism/sample_scores.tsv",
        "diagnostic": args.data_dir/"mechanism/within_method_diagnostics.tsv",
        "robust": args.data_dir/"robustness/performance_metrics.tsv",
        "robust_scores": args.data_dir/"robustness/sample_scores.tsv",
    }
    data = {name: read(path) for name, path in source_files.items()}
    fragments = fragment_data(args.fragment_data, data["robust_scores"])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    jobs = [
        ("Experiment1_Detection_performance", lambda: experiment1(data["primary"])),
        ("Experiment2_Four_method_distributions", lambda: experiment2(data["mechanism_scores"])),
        ("Experiment3_Depth_purity_direction", lambda: experiment3(data["robust"])),
    ]
    if not args.no_heatmap:
        jobs.append(("Supplementary_Depth_purity_heatmap", lambda: supplementary_heatmap(data["robust"])))
    for name, make in jobs:
        fig = make()
        for suffix in dict.fromkeys(args.formats):
            path = args.out_dir/f"{name}.{suffix}"
            fig.savefig(path, dpi=args.dpi, bbox_inches="tight", pad_inches=.16)
            print(path)
        plt.close(fig)
    tables(args.out_dir, data["primary"], data["mechanism"], data["diagnostic"], data["robust"], fragments)
    (args.out_dir/"Figure_panel_definitions.md").write_text(PANEL_GUIDE, encoding="utf-8")
    # Record data identity without modifying the underlying TSV files.
    sources = list(source_files.values())+[args.fragment_data]
    (args.out_dir/"Source_checksums.txt").write_text("\n".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.resolve()}" for path in sources)+"\n", encoding="utf-8")
    print("Three Markdown results tables and panel definitions written.")


if __name__ == "__main__":
    main()
