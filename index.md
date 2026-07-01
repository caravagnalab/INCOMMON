# INCOMMON

The `INCOMMON` package implements the statistical model described in
[Calonaci et al; medRxiv
(2024)](https://doi.org/10.1101/2024.05.13.24307238) to infer the copy
number and multiplicity of somatic mutations from tumor-only read count
data. Mutations are classified as either Tier-1 (present in 100% of
cells) without copy-number alterations (heterozygous mutant diploid),
with loss of heterozygosity, copy-neutral loss of heterozygosity,
amplification, or Tier-2 (subclonal or with high ploidy and low
multiplicity), providing a fast genome-interpretation framework for
large-scale tumor-only sequencing cohorts.

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

If you want to run INCOMMON classification
([`classify()`](caravagnalab.github.io/INCOMMON/reference/classify.md))
yourself, you also need [`cmdstanr`](https://mc-stan.org/cmdstanr/) and
a working Stan toolchain, which are not installed automatically:

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
