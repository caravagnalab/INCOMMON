# INCOMMON

``` r

library(INCOMMON)
#> Warning: replacing previous import 'cli::num_ansi_colors' by
#> 'crayon::num_ansi_colors' when loading 'INCOMMON'
```

INCOMMON is a tool for the INference of COpy number and Mutation
Multiplicity in ONcology. INCOMMON infers the copy number and
multiplicity of somatic mutations from tumor-only read count data, and
can be applied to tumour-only samples in an efficient and fast way. For
each mutation in a sample, INCOMMON computes a probability distribution
over all the possible combinations of total copy number and
multiplicity, and choses the one with maximum probability.

In addition, INCOMMON offers a genome interpretation framework, in which
the tumour genome is classified based on the mutant dosage of oncogenes
and tumour suppressor genes. Clinical outcome analysis (e.g. survival
analysis and organotropism) based on such complex genotypes can be
performed with functions integrated in the package.

> The INCOMMON model is designed to work with high coverage sequencing
> data such as targeted panels but, in principle, it can be used with
> any sequencing assay. INCOMMON is helpful also to analyse sequencing
> data from tumor-only assays, in paricular when alignment files (fastq,
> sam/bam, etc.) are not availble. However, if one can access
> higher-resolution whole-exome or whole-genome assays, specific
> [deconvolution methodologies should be
> used](https://caravagnalab.github.io/mobster,%20https://github.com/caravagnalab/CNAqc).

## The INCOMMON mutation copy-number caller

INCOMMON is a Bayesian model that can infer, for a tumour mutation,
$`(i)`$ the total copy number at the mutant locus, and $`(ii)`$ the
mutation multiplicity (i.e., the number of DNA copies that harbour the
mutation). This information provides the allele-specific configuration
of the mutant locus.

INCOMMON takes input read counts data for $`n`$ mutations
$`X=\{x_1, \ldots, x_n\}`$ to develop the joint likelihood}
``` math
\begin{equation}
{p(X\mid \Theta) =  \prod_{x_i \in X} p(x_i \mid \Theta)
=
\prod_{x_i\in X} p(d_i \mid \Theta)\, p(r_i \mid d_i, \Theta)
\,.}
\end{equation}
```

For every mutation $`x_i=\langle r_i, d_i\rangle`$, INCOMMON uses the
number of reads $`r_i`$ with the alternative allele and the total reads
$`d_i`$ (depth of sequencing). The model infers two sample-specific
parameters, ($`i`$) the sample purity ($`0 < \pi \leq 1`$), and ($`ii`$)
the rate of reads per chromosome copy ($`\eta \in \mathbb{R}^+`$), and
two mutation-specific parameters, ($`iii`$) the tumour total copy number
at the locus ($`k_i\in \mathbb{Z}^+`$) and ($`iv`$) the mutation
multiplicty ($`m_i\in \mathbb{Z}^+`$, $`m_i\leq k_i`$).

The model adopted by INCOMMON relates copy number to coverage linearly,
so the expected number of reads for $`k`$ chromosome copies is
$`k\eta`$.

INCOMMON links the number of reads to the multiplicity and total copy
number considering tumour/ normal admixing
``` math
\begin{align}
p(d_i\mid \pi, \eta, k_i) &
= \text{Poisson}(d_i \mid \lambda)
&& 
\lambda  = (1-\pi)2\eta + \pi \eta k_i  \label{eq:likelihood_pois}\\
p(r_i\mid  \pi, k_i, m_i) &
= \text{Binomial}(r_i \mid d_i, \varphi)
&&
\varphi = \frac{m_i\pi}{2\left(1-\pi\right)+k_i\pi}\label{eq:likelihood_binom}
\end{align}
```

The sequencing depth follows a Poisson distribution with an expected
number of reads $`\lambda`$ defined by combining tumour and normal
readouts. Given the depth, the number of mutant reads follows a Binomial
distribution with success probability $`\varphi`$ determined by the
mixing of tumour and normal success rates.

INCOMMON uses Markov Chain Monte Carlo (implemented in `stan`) sampling
to estimate a posterior distribution over $`m`$, $`k`$, $`\eta`$ and
$`\pi`$. To leverage the massive amount of public WGS data of human
tumours and gain precision with targeted assays, by default INCOMMON
uses [a biologically informed prior distribution for copy number and
multiplicity configurations from the PCAWG and HMF
cohorts](caravagnalab.github.io/INCOMMON/articles/a1_priors_k_m.md).

To support orthogonal estimation of tumour purity (e.g. from
histopathological evaluation) but resist potential error in the input
estimate, INCOMMON centres a prior around the purity measurements
provided with each sample. After posterior inference, INCOMMON uses
posterior predictive checks to monitor the discrepancy between observed
values and inferred posterior distributions on $`k`$, $`m`$, $`\eta`$
and $`\pi`$.

## Input format

The input required for INCOMMMON classification consists of two data
tables:

- `genomic_data`: a data table of annotated mutations with columns
  indicating, for each mutation, the sample name `sample`, mutant
  chromosome `chr`, start position `from`, end position `to`, reference
  allele `ref`, alternative allele `alt`, sequencing depth `DP`, number
  of reads with variant `NV`, variant allele frequency `VAF`, mutant
  gene name `gene` Hugo Symbol, and possibly the protein sequence of the
  variant in HGVS recommended format (preferably 1-letter amino-acid
  code `HGVSp_Short`).

- `clinical_data`: a data table of clinical data with matched sample
  names `sample` and purity `purity` (required), and clinical features
  like tumor type (ONCOTREE code) `tumor_type` (required for tumor
  specific priors), survival data such as `OS_STATUS` and time
  `OS_MONTHS` (required for survival analysis), metastasis data such as
  `SAMPLE_TYPE` (Primary or Metastasis), number of metastases
  `MET_COUNT` (required for metastatic propensity analysis) and
  metastatic site `METASTATIC_SITE` (required for metastatic tropism
  analysis), plus any other useful covariate.

- `gene_roles`: a data table reporting gene names `gene` and associated
  roles `gene_role` (“oncogene” or “TSG”). INCOMMON provides a set of
  gene roles extracted from the COSMIC Cancer Gene Census (v98) as
  default.

The input for downstream analysis is checked and cast in the expected
format through the function `init`.

INCOMMON provides data from the publicly available MSK-MetTropism cohort
in the correct format. The following example shows how this input is
pre-processed by INCOMMON:

``` r

data(MSK_genomic_data)
data(MSK_clinical_data)
data(cancer_gene_census)

x = init(
  genomic_data = MSK_genomic_data,
  clinical_data = MSK_clinical_data,
  gene_roles = cancer_gene_census
)
#> ── INCOMMON - Inference of copy number and mutation multiplicity in oncology ───
#> 
#> ── Genomic data ──
#> 
#> ✔ Found 25659 samples, with 224939 mutations in 491 genes
#> ! No read counts found for 1393 mutations in 1393 samples
#> ! Gene name not provided for 1393 mutations
#> ! 201 genes could not be assigned a role (TSG or oncogene)
#> 
#> ── Clinical data ──
#> 
#> ℹ Provided clinical features:
#> ✔ sample (required for classification)
#> ✔ purity (required for classification)
#> ✔ tumor_type
#> ✔ OS_MONTHS
#> ✔ OS_STATUS
#> ✔ SAMPLE_TYPE
#> ✔ MET_COUNT
#> ✔ METASTATIC_SITE
#> ✔ MET_SITE_COUNT
#> ✔ PRIMARY_SITE
#> ✔ SUBTYPE_ABBREVIATION
#> ✔ GENE_PANEL
#> ✔ SEX
#> ✔ TMB_NONSYNONYMOUS
#> ✔ FGA
#> ✔ AGE_AT_SEQUENCING
#> ✔ RACE
#> ✔ Found 25257 matching samples
#> ✖ Found 513 unmatched samples

print(x)
#> ── [ INCOMMON ]  223546 PASS mutations across 24266 samples,
#> with 490 mutant gen
#> ℹ Average sample purity: 0.4
#> ℹ Average sequencing depth: 660
#> # A tibble: 223,546 × 27
#>    sample    tumor_type purity chr     from     to ref   alt      DP    NV   VAF
#>    <chr>     <chr>       <dbl> <chr>  <dbl>  <dbl> <chr> <chr> <int> <int> <dbl>
#>  1 P-0028912 CHOL          0.3 chr17 7.58e6 7.58e6 G     A       837   133 0.159
#>  2 P-0028912 CHOL          0.3 chr6  1.12e8 1.12e8 -     A       698   141 0.202
#>  3 P-0028912 CHOL          0.3 chrX  5.32e7 5.32e7 G     A       832    85 0.102
#>  4 P-0003698 BLCA          0.2 chr17 7.58e6 7.58e6 C     A       437   109 0.249
#>  5 P-0003698 BLCA          0.2 chr3  4.99e7 4.99e7 C     A       591    86 0.146
#>  6 P-0003698 BLCA          0.2 chr5  1.49e8 1.49e8 C     T       360    36 0.1  
#>  7 P-0003698 BLCA          0.2 chr13 3.29e7 3.29e7 G     C      1027   162 0.158
#>  8 P-0003698 BLCA          0.2 chr13 3.29e7 3.29e7 G     C      1021   182 0.178
#>  9 P-0003698 BLCA          0.2 chr19 1.11e7 1.11e7 G     T       573    98 0.171
#> 10 P-0003698 BLCA          0.2 chr22 4.15e7 4.15e7 G     A       416    45 0.108
#> # ℹ 223,536 more rows
#> # ℹ 16 more variables: gene <chr>, gene_role <chr>, OS_MONTHS <dbl>,
#> #   OS_STATUS <dbl>, SAMPLE_TYPE <chr>, MET_COUNT <dbl>, METASTATIC_SITE <chr>,
#> #   MET_SITE_COUNT <dbl>, PRIMARY_SITE <chr>, SUBTYPE_ABBREVIATION <chr>,
#> #   GENE_PANEL <chr>, SEX <chr>, TMB_NONSYNONYMOUS <dbl>, FGA <dbl>,
#> #   AGE_AT_SEQUENCING <dbl>, RACE <chr>
```

## Gene mutant dosage

Downstream of INCOMMON classification, the tumour genome can be
interpreted in terms of the [gene mutant
dosage](caravagnalab.github.io/INCOMMON/articles/a3_mutant_dosage.md) of
tumor suppressor genes (typically through mutations with LOH) and
enhanced activation of oncogenes (typically through mutations with
mutant copy gain).

Mutant dosage classes (low, balanced and high) are derived from the
fraction of alleles carrying the mutation (FAM)

``` math
\begin{equation}\label{e:expected_fma}
    \text{FAM}_{g,t}=\sum\limits_{k=1}^{k_{max}}\sum_{m=1}^k\frac{m}{k}p(m,k\:|\: X_{g,t})\, ,
\end{equation}
```

which is derived from the full posterior distribution of
$`p(m,k\:|\: X_{g,t})`$, computed by INCOMMON.

## Survival analysis

If patients’ survival status and time are provided as features in the
clinical data table `clinical_table`, [survival
analysis](caravagnalab.github.io/INCOMMON/articles/a4_survival_analysis.md)
can be performed. Downstream of classification, INCOMMON can stratify
patients based on the mutant dosage of a TSG or an oncogene of interest.

INCOMMON provides the following functions, dedicated to fitting survival
data:

- `kaplan_meier_fit` uses the Kaplan-Meier estimator to fit survival
  data from patients stratified with respect to the status of a specific
  `tumor_type` and `gene` (low, balanced or high mutant dosage)
- `cox_fit` uses a Cox proportional hazard ratio model to fit survival
  data. In addition to arguments `tumor_type` and `gene`, it accepts
  other `covariates`, given they are provided in the `clinical_table`.

## Metastatic patterns

If information about metastatisation is provided, such as type of the
sample (primary tumor or metastasis), whether patients are metastatic or
not, and sites of metastatisation for primary tumors, in the clinical
data table `clinical_table`, [analysis of metastatic propensity and
tropism](caravagnalab.github.io/INCOMMON/articles/a5_metastasis_analysis.md)
based on INCOMMON classification and genome interpretation can be
performed.

INCOMMON provides the following functions for analysis of metastases:

- `met_propensity`: uses a logistic regression test to compute the odds
  ratio (OR) of metastatisation between patients identified by the
  mutant dosage of a `gene` (low, balanced or high) for specific types
  of primary tumors `tumor_type`.

- `met_tropism`: uses a logistic regression test to compute the odds
  ratio (OR) to metastatise to a specific site `METASTATIC_SITE` between
  patients identified by the two mutational statuses of a `gene` (Mutant
  TSG with versus without LOH and Mutant oncogene with versus without
  amplification) for specific types of primary tumors `tumor_type`.
