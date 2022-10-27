purity_from_fit = function(n_binomials, peaks, purity, eps) {
  out = NA
  if (n_binomials == 3)
    out = peaks[2] * 2
  if (n_binomials == 1)
    out = peaks[1] * 2
  if (n_binomials == 2 & !is.na(purity) & min(abs(peaks[2] * 2 - purity), abs(peaks[1] * 2 - purity)) <= eps)
    out = peaks[which.min(c(abs(peaks[2] * 2 - purity), abs(peaks[1] * 2 - purity)))] * 2
  if (is.na(out)) {
    cli::cli_alert("Purity could not be estimated reliably. Check results of the fit in the output figure")
  }
  return(out)
}