purity_from_fit = function(n_binomials, peaks, purity, eps){
  out = case_when(
    n_binomials == 3 ~ peaks[2] * 2,
    n_binomials == 1 ~ peaks[1] * 2,
    n_binomials == 2 &
      min(abs(peaks[2] * 2 - purity), abs(peaks[1] * 2 - purity)) <= eps ~ peaks[which.min(c(abs(peaks[2] *
                                                                                                   2 - purity), abs(peaks[1] * 2 - purity)))] * 2,
    TRUE ~ NA %>% as.double()
  )
  if(is.na(out)){
    cli::cli_alert("Purity could not be estimated reliably. Check results of the fit in the output figure")
  }
  return(out)
}