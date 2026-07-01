# 6. Interactive analysis with the INCOMMON web app

All the analyses described in the other vignettes
([classification](caravagnalab.github.io/INCOMMON/articles/a2_classify_mutations.md),
[gene mutant
dosage](caravagnalab.github.io/INCOMMON/articles/a3_mutant_dosage.md),
[survival
analysis](caravagnalab.github.io/INCOMMON/articles/a4_survival_analysis.md)
and [metastatic pattern
analysis](caravagnalab.github.io/INCOMMON/articles/a5_metastasis_analysis.md))
are also available without writing any R code, through the INCOMMON web
app:

**<https://ncalonaci.shinyapps.io/incommon/>**

The app runs the same underlying INCOMMON model used by this package,
and is meant for users who want to explore the method interactively or
classify a dataset without setting up an R/`cmdstanr` environment.

With the app you can, briefly:

- upload your own `genomic_data` and `clinical_data` tables (in the same
  format expected by
  [`init()`](caravagnalab.github.io/INCOMMON/reference/init.md), see the
  [Get started](caravagnalab.github.io/INCOMMON/articles/INCOMMON.md)
  vignette) or explore the bundled MSK-MetTropism example cohort;
- run INCOMMON classification of mutation copy number and multiplicity
  directly from the browser, without installing `cmdstanr`/Stan locally;
- browse the resulting genome interpretation (gene mutant dosage
  classes) and the diagnostic/posterior plots produced by the model;
- interactively perform survival and metastatic-tropism analyses on the
  classified cohort;
- download the classification results and plots for use outside the app.

For full control over the analysis (custom priors, scripting,
integration into a larger pipeline), use the R package directly as
described in the other vignettes; the app is best suited to quick,
exploratory use.
