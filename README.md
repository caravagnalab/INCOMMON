
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/INCOMMON/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/caravagnalab/INCOMMON/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/caravagnalab/INCOMMON/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/caravagnalab/INCOMMON/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Citations](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fapi.semanticscholar.org%2Fgraph%2Fv1%2Fpaper%2FDOI%3A10.1101%2F2024.05.13.24307238%3Ffields%3DcitationCount&query=%24.citationCount&label=citations&color=blue)](https://www.semanticscholar.org/paper/2d5318baf1a79d491822e5d98a1e06028f521fbc)
[![GitHub
stars](https://img.shields.io/github/stars/caravagnalab/INCOMMON?style=flat&color=yellow)](https://github.com/caravagnalab/INCOMMON/stargazers)
<!-- badges: end -->

# INCOMMON <a href='https://caravagnalab.github.io/INCOMMON'><img src='man/figures/logo.png' align="right" height="139" /></a>

The `INCOMMON` package implements the statistical framework described in
[Calonaci et al; medRxiv
(2024)](https://doi.org/10.1101/2024.05.13.24307238) to infer mutation
copy number and multiplicity directly from tumor-only clinical targeted
sequencing data, without requiring matched normal samples. This
information is used to derive gene mutant dosage, an emergent property
of the interplay between somatic mutations and copy-number alterations
that is overlooked by standard binary mutant/wild-type models. INCOMMON
stratifies patients by mutant dosage of oncogenes and tumor suppressor
genes to identify biomarkers of prognosis and metastatic tropism.

INCOMMON is also available as a user-friendly
[ShinyApp](https://ncalonaci.shinyapps.io/incommon/) (see [this
vignette](https://caravagnalab.github.io/INCOMMON/articles/a6_shiny_app.html)
for a brief overview of what it offers).

You can download the results of our analysis from
[Zenodo](https://zenodo.org/records/12547426).

#### Citation

[![](https://img.shields.io/badge/doi-10.1101/2024.05.13.24307238-red.svg)](https://doi.org/10.1101/2024.05.13.24307238)

If you use `INCOMMON`, please cite our preprint:

- Calonaci, N., Krasniqi, E., Čolić, D. et al. *Gene mutant dosage is
  associated with prognosis and metastatic tropism in 60,000 clinical
  cancer samples.* medRxiv (2024).

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/INCOMMON/-yellow.svg)](https://caravagnalab.github.io/INCOMMON)

------------------------------------------------------------------------

### Installation

You can install INCOMMON from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/INCOMMON")
```

If you want to run INCOMMON classification (`classify()`) yourself, you
also need [`cmdstanr`](https://mc-stan.org/cmdstanr/) and a working Stan
toolchain, which are not installed automatically:

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
```

Alternatively, the [INCOMMON web
app](https://ncalonaci.shinyapps.io/incommon/) (see [this
vignette](https://caravagnalab.github.io/INCOMMON/articles/a6_shiny_app.html))
runs the same model in your browser, with no local Stan installation
required.

------------------------------------------------------------------------

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory. University of Trieste, Italy.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
