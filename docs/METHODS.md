# Reproducible methods

MSIpath ranks normal reference donors using ancestry-marker and flanking-SNP likelihoods, retains 20 donors (up to 40 haplotypes), then infers local copying paths and marginal genotype weights.

At each locus, the null hypothesis uses a normal length model weighted by genotype and allelic balance. The alternative hypothesis changes one allele by a fixed length shift shared across all observations at that locus, jointly marginalizing over genotype, allelic balance, the altered allele, shift size, and altered-molecule fraction once per locus. Expansion and contraction directions receive equal prior probability. The fixed shift-size distribution spans 1–40 repeat units and mixes a short-range geometric distribution with weight 0.7 and a broad-range geometric distribution with weight 0.3.

The altered-molecule fraction prior is a truncated Beta(1,9) distribution integrated using 12-point quadrature, with an upper bound constrained by both the selected allele's proportion and 0.5. The observation model retains mapping uncertainty. Only fragments that permit exact determination of STR length enter the likelihood used for scoring; each eligible locus requires at least 3 such exact fragments and at least 5 fragments in total.

The sample score is the sum of locus log Bayes factors individually clipped to [0,10], requiring at least 5 eligible loci. Ranking the query score against independent normal sample scores gives:

$$
p_{\mathrm{calibrated}} = \frac{1 + \text{number of normal scores greater than or equal to the query score}}{N+1}
$$

Here, $N$ is the number of normal samples used for calibration. The fixed significance level $\alpha$ is 0.05. A query sample with insufficient evidence is labeled `UNCALLABLE`; a scored query sample without a calibration bank is labeled `UNCALIBRATED`.

The bundled normal error profiles were fitted from 40 reference-only donors using maximum a posteriori (MAP) normal-genotype residual categories and pooled shrinkage.

The study simulates PCR-free PE150 paired-end reads from target regions on chr19 and chr22, using 200 selected STR loci, phased SNP/STR templates, ART sequencing errors, and a primary sequencing depth of 30×. The robustness study includes depths of 10×, 20×, 30×, and 60×. Tumor purity $p$ is defined under clonal heterozygous alterations in diploid, copy-neutral cells, so the altered-molecule fraction equals $p/2$. The parameter $b$ is the fraction of altered loci. The remaining loci can stay stable even in a tumor sample; $b$ and $p$ are distinct parameters.

Calibration and performance summaries use independent normal cohorts sequenced at the same depth as the query samples. Native tool calls and results obtained using the frozen 30× calibration scheme are retained as secondary columns. AUROC (area under the receiver operating characteristic curve) and AUPRC (area under the precision–recall curve) compare each positive condition against ordinary normal samples; challenging normal samples have separate false-positive results. Confidence intervals use donor-level resampling or Wilson binomial intervals, as specified in the tables.

Resource measurements include fresh process startup, prebuilt model loading, BAM evidence extraction, inference, calibration, and output. They exclude alignment shared across methods, simulation, one-time model or baseline training, and output archiving. Each method processes the same 10 tumor- or normal-state samples from 5 donors, with 3 fresh-process repeats and affinity to a single CPU. MANTIS2 additionally requires matched normal BAM files. All five methods share 15 supported loci; a separate 200-locus comparison omits MSIsensor2.

Agreement in scores and callable loci was checked for all 270 executions. The published resource evaluation has a small sample size and retains uncallable executions.
