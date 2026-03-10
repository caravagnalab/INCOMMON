# Package index

## Available data

Datasets included in the package.

- [`MSK_genomic_data`](caravagnalab.github.io/INCOMMON/reference/MSK_genomic_data.md)
  : Genomic data of the MSK-MetTropism cohort
- [`MSK_clinical_data`](caravagnalab.github.io/INCOMMON/reference/MSK_clinical_data.md)
  : Genomic data of the MSK-MetTropism cohort
- [`priors_pcawg_hmf`](caravagnalab.github.io/INCOMMON/reference/priors_pcawg_hmf.md)
  : Priors over (k,m) configurations from PCAWG and HMF
- [`priors_eta`](caravagnalab.github.io/INCOMMON/reference/priors_eta.md)
  : Priors over the rate of reads per chromosome copy from MSK-MET
- [`cancer_gene_census`](caravagnalab.github.io/INCOMMON/reference/cancer_gene_census.md)
  : Gene roles from COSMIC Cancer Gene Census
- [`MSK_PAAD_output`](caravagnalab.github.io/INCOMMON/reference/MSK_PAAD_output.md)
  : Data from the MSK-MetTropism cohort classified with INCOMMON

## INCOMMON classification

Main functions to classify mutations

- [`compute_eta_prior()`](caravagnalab.github.io/INCOMMON/reference/compute_eta_prior.md)
  : Evaluate a prior distribution over the rate of reads per chromosom
  copy from the data.

- [`priors_k_m()`](caravagnalab.github.io/INCOMMON/reference/priors_k_m.md)
  :

  Getter for class `'INCOMMON'`.

- [`subset_sample()`](caravagnalab.github.io/INCOMMON/reference/subset_sample.md)
  : Subset an INCOMMON object by sample ID.

- [`init()`](caravagnalab.github.io/INCOMMON/reference/init.md) :

  Prepare input for analyses with `'INCOMMON'`.

- [`classify()`](caravagnalab.github.io/INCOMMON/reference/classify.md)
  : Classify mutations using a Beta-Binomial model-based test.

## Visualisation

Functions to plot INCOMMON classification results

- [`print(`*`<INCOMMON>`*`)`](caravagnalab.github.io/INCOMMON/reference/print.INCOMMON.md)
  :

  Print for class `'INCOMMON'`.

- [`plot_prior()`](caravagnalab.github.io/INCOMMON/reference/plot_prior.md)
  : Visualize prior distribution for a gene (tumor-specific or
  pancancer).

- [`plot_eta_prior()`](caravagnalab.github.io/INCOMMON/reference/plot_eta_prior.md)
  : Visualise the prior distribution over the rate of reads per
  chromosome copy.

- [`plot_purity_prior()`](caravagnalab.github.io/INCOMMON/reference/plot_purity_prior.md)
  : Visualise the prior distribution on sample purity

- [`plot_posterior_k_m()`](caravagnalab.github.io/INCOMMON/reference/plot_posterior_k_m.md)
  : Visualise the posterior distribution on (k,m) configurations.

- [`plot_class_fraction()`](caravagnalab.github.io/INCOMMON/reference/plot_class_fraction.md)
  : Visualize frequency distribution of INCOMMON classes.

- [`plot_survival_analysis()`](caravagnalab.github.io/INCOMMON/reference/plot_survival_analysis.md)
  : Visualise survival analysis based on INCOMMON classes.

- [`plot_met_volcano()`](caravagnalab.github.io/INCOMMON/reference/plot_met_volcano.md)
  : Visualize metastatic propnesity odds ratio in a volcano plot
  fashion.

- [`plot_tropism()`](caravagnalab.github.io/INCOMMON/reference/plot_tropism.md)
  : Visualize metastatic propnesity odds ratio in a volcano plot
  fashion.

## Getter Functions

Functions to extract model information

- [`parameters()`](caravagnalab.github.io/INCOMMON/reference/parameters.md)
  :

  Getter for class `'INCOMMON'`.

- [`show_FAM()`](caravagnalab.github.io/INCOMMON/reference/show_FAM.md)
  : Get the fraction of alleles with the mutation (FAM) values for a
  gene and cancer type.

- [`posterior_k_m()`](caravagnalab.github.io/INCOMMON/reference/posterior_k_m.md)
  : Get the posterior distribution over (k,m) configurations for a
  specific mutation

## INCOMMON Model

Functions to compute model likelihood and posterior distributions

- [`mutant_dosage_classification()`](caravagnalab.github.io/INCOMMON/reference/mutant_dosage_classification.md)
  : Group patients by gene mutant mutant dosage using gene-role specific
  thresholds.
- [`compute_expectations()`](caravagnalab.github.io/INCOMMON/reference/compute_expectations.md)
  : Compute the expectation value of k, m and FAM from the model full
  posterior distribution.

## Survival analysis with INCOMMON classes

Functions to perform survival analysis with patient stratificaiton based
on INCOMMON classes

- [`kaplan_meier_fit()`](caravagnalab.github.io/INCOMMON/reference/kaplan_meier_fit.md)
  : Fit Kaplan-Meier survival model stratifying patients based on
  INCOMMON classes, for a single gene and tumour type.
- [`cox_fit()`](caravagnalab.github.io/INCOMMON/reference/cox_fit.md) :
  Fit multivariate Cox regression model based on INCOMMON classes.
- [`plot_survival_analysis()`](caravagnalab.github.io/INCOMMON/reference/plot_survival_analysis.md)
  : Visualise survival analysis based on INCOMMON classes.

## Metastatic pattern analysis with INCOMMON classes

Functions to perform metastatic propensity and tropism analysis with
INCOMMON interpreted tumor genomes

- [`met_propensity()`](caravagnalab.github.io/INCOMMON/reference/met_propensity.md)
  : Logistic regression of metastatic propensity based on INCOMMON
  classes.
- [`met_tropism()`](caravagnalab.github.io/INCOMMON/reference/met_tropism.md)
  : Logistic regression of metastatic tropism based on INCOMMON classes.
