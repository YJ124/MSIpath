# MSIpath

Repository: https://github.com/YJ124/MSIpath

MSIpath detects microsatellite instability using an SNP-recommended personalized reference, local haplotype-path inference, and a marginal single-altered-allele likelihood. The released `msipath` command runs this model directly. The sample score sums positive log Bayes factors clipped at 10 and is calibrated against independent normal donors.

## Install

Python 3.10+ on Linux or Windows via WSL is recommended. Install `samtools` for BAM/CRAM input; SAM and saved-evidence examples do not require it.

```bash
python -m venv .venv
source .venv/bin/activate
pip install '.[dev]'
pytest
msipath --version
```

## Run from an alignment

```bash
msipath call --bam sample.bam --reference GRCh38.fa \
  --panel examples/research_panel --normal-profile examples/normal_profiles.json \
  --calibration examples/calibration_panel200.json \
  --sample-id sample01 --donor-id donor01 --output sample01-results
```

Inputs must use the same assembly and chromosome names as the panel. The included panel, normal profiles and calibration banks reproduce the chr19/chr22, 30x, PCR-free simulation study. They are not universal clinical assets. Use independently prepared normal calibration data for a new sequencing protocol or locus scope.

## Self-contained example

```bash
python examples/run_demo.py --output demo-output
```

This example uses saved SNP counts and STR observations; it recomputes recommendation, local paths and the released likelihood. It excludes BAM extraction and should not be used as a BAM runtime benchmark.

## Independent normal calibration

Run `msipath call` on each calibration normal without `--calibration`. It returns a score and `UNCALIBRATED`, never an invented MSI/MSS call. Prepare a tab-separated manifest with `donor`, `role` (all `calibration`), `label` (all `0`), and `result` (path to each result.json).

```bash
msipath calibrate --manifest calibration_normals.tsv --output calibration.json
```

At least 20 distinct usable normal donors are required; the paper uses a 40-donor cohort. Reference/calibration/query overlap and mismatched model/profile/panel/locus signatures are rejected. Calibration samples must also be independent of normal-profile training, as documented in the study protocol.

## Reference building

`msipath build-reference --config configs/reference.example.yaml --output new-panel` builds a schema-v2 panel from jointly phased SNP/STR sequence data, locus registries and graph pointers. It does not assume that independently phased VCFs are linked merely because sample IDs match. See the configuration example and `docs/METHODS.md`.

## Outputs

`result.json` contains `score`, `callable_loci`, `calibrated_p`, and `call` (MSI, MSS, UNCALLABLE or UNCALIBRATED). Locus log Bayes factors, recommendation, observations and GFA graphs are retained. Locus Bayes factors are not posterior MSI probabilities or multiple-testing-adjusted p-values.

## Paper artifacts

`paper/` contains the final five-method tables and figures, a standalone HTML report and a printable PDF. `analysis/build_report.py` rebuilds these from the included tables. `analysis/export_results.py --workspace /path/to/original/experiments` imports the original frozen results without changing them. Only MSIpath, MSIsensor-pro, MSIsensor2, mSINGS and MANTIS2 appear in the method comparisons.

Tumor purity is the public simulation axis. For this study's clonal heterozygous, diploid copy-neutral model, the altered-molecule fraction is half the tumor purity. Conversion changes the numerical values; observed scores, calls and confidence intervals are not modified. Burden is the fraction of altered STR loci.

## Validation and scope

The release is checked against frozen full-study outputs, including a fresh BAM run. See `validation/`. The study reuses development-era simulation donors; it is not blind clinical validation. Difficult inherited STR alleles remain a major false-positive failure mode, and MSIpath uses more time/memory than the four comparator tools in the measured small-panel benchmark. All these results remain in the final report.

External comparator binaries, raw sequencing data and genome FASTA files are not bundled. Citation/source links are in `docs/SOURCES.md`; public 1000 Genomes donor identifiers in the research panel are retained for provenance. This repository contains the released MSIpath implementation and simulation-study artifacts.

## Complete end-to-end purity grid

`paper/data/end_to_end_purity_grid/` contains the completed 20%, 40%, 60%, and 80% purity evaluation. Each burden/purity condition uses the same 200 test donors; repeated conditions are not independent patients. The original `paper/data/end_to_end/` and summary report are retained as the earlier experiment snapshot and do not include all later purity supplements.

The updated figure is in `paper/figures/end_to_end_purity_grid/`. To reproduce it:

```bash
pip install '.[report]'
python analysis/plot_endtoend_purity_arial.py --font-dir /path/to/Arial/fonts
```

Arial fonts must be installed separately and are not redistributed here. Every bar is checked against the included individual sample calls. Unresolved samples remain in the denominator; error bars are Wilson 95% confidence intervals. Panels A/B use 200 loci with burdens 10%/30%; panels C/D use 15 common loci with burdens 10%/30%.

No patient WBC files or patient-level WBC results are included. This repository contains simulated observations and public-reference donor identifiers.
