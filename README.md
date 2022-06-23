
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TAPACLOTH <a href="caravagnalab.github.io/TAPACLOTH"><img src="man/figures/logo.png" align="right" height="139" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/caravagnalab/TAPACLOTH/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/caravagnalab/TAPACLOTH/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

TAPACLOTH is a tool that uses targeted sequencing data to classify
somatic mutations as either subclonal, clonal heterozygous or clonal
with loss of heterozygosity. Two types of tests are available: an exact
one based on Binomial or Beta-Binomial modelling of sequencing read
counts, and approximate based on allelic frequencies. The former is
parametrised by a p-value. Additionally, TAPACLOTH can infer sample
purity, an information used to carry out variant classification. Since
the two tests are independent, one can first estimate sample purity and
then use it as input for the classification task.

#### Help and support

## [![](https://img.shields.io/badge/GitHub%20Pages-https://caravagnalab.github.io/TAPACLOTH/-yellow.svg)](https://caravagnalab.github.io/TAPACLOTH)

## Installation

You can install the development version of TAPACLOTH from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("caravagnalab/TAPACLOTH")
```

#### Copyright and contacts

Cancer Data Science (CDS) Laboratory.

[![](https://img.shields.io/badge/CDS%20Lab%20Github-caravagnalab-seagreen.svg)](https://github.com/caravagnalab)
[![](https://img.shields.io/badge/CDS%20Lab%20webpage-https://www.caravagnalab.org/-red.svg)](https://www.caravagnalab.org/)
