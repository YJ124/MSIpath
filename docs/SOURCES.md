# Sources and external software

- ART: Huang et al., *ART: a next-generation sequencing read simulator*, Bioinformatics (2012), https://doi.org/10.1093/bioinformatics/btr708.
- Population STR reference: EnsembleTR population resource, https://pmc.ncbi.nlm.nih.gov/articles/PMC10593948/.
- MSIsensor-pro: https://github.com/xjtu-omics/msisensor-pro.
- MSIsensor2: https://github.com/niu-lab/msisensor2.
- mSINGS: this study used the checkout at https://github.com/sheenamt/msings; its exact commit is recorded in `validation/environment.json`.
- MANTIS2: https://github.com/nh13/MANTIS2. The exact upstream checkout used in the study is recorded in `validation/environment.json`; the study uses paired-normal mode.

External software and genomic reference distributions are not vendored into this package. Dependency and upstream licenses apply to their own distributions. Public reference donor IDs are preserved in the included research panel for exclusion and reproducibility.

See `paper/data/provenance.json` for hashes of the frozen result tables. Only display names and simulation purity units were transformed. Original private workspace results were kept unchanged.
