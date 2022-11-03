purity_from_fit = function(n_binomials, peaks, purity, eps) {
  out = NA
  if (n_binomials == 3)
    out = peaks[2] * 2
  if (n_binomials == 1)
    out = peaks[1] * 2
  if (n_binomials == 2){
    purity_hyp = c(min(peaks[1] * 2, 1), min(peaks[2] * 2, 1)) # in this case we have two equivalent putative purities 
    if(!is.na(purity) & min(abs(purity_hyp-purity)) <= eps) # check if at least one is close to the input one
    out = purity_hyp[which.min(abs(purity_hyp-purity))] # choose the estimate that is closest to input 
  if (is.na(out)) {
    cli::cli_alert("Purity could not be estimated reliably. Check results of the fit in the output figure")
  }
  return(out)
}