# INCOMMON 0.1.0

* Package-hygiene pass, mirroring the recent revolver/ctree/VIBER cleanup:
  - Removed the redundant `Author:` field from `DESCRIPTION` (superseded by `Authors@R`).
  - Fixed duplicated `cre` role in `Authors@R` (Nicola Calonaci is the maintainer; Giulio Caravagna is `aut`), fixed a typo'd email address, and added the missing ORCID.
  - Rewrote `Title` to follow the `PKG: Title Case` convention.
  - Fixed `URL` to use full `https://` links and include the GitHub repository.
  - Guarded the unconditional `cmdstanr::cmdstan_model()` call in `get_stan_model()` with `requireNamespace()`, since `cmdstanr` is a `Suggests`-only dependency (not on CRAN).
  - Modernized GitHub Actions workflows: `actions/checkout@v4`, concurrency groups to cancel superseded runs, a simplified release-only ubuntu+macos `R-CMD-check` matrix, and a `pkgdown` workflow that deploys via `pkgdown::deploy_to_branch()` instead of a third-party action (cmdstan installation steps preserved in both).
  - Untracked stray files that had been committed by mistake (`.DS_Store`, `.Rhistory`, an empty stray file) and renamed the leftover `TAPACLOTH.Rproj` to `INCOMMON.Rproj`.
  - Resized the package logo to hex-sticker proportions (240x282).
  - Added `pkgdown`/`Lifecycle` badges and a preprint DOI badge to `README`.
  - Bumped version to 0.1.0.
