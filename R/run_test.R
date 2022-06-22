run_test = function(nv, null_model){
  if(nv > null_model$nv[2]) return("Clonal LOH")
  if(nv > null_model$nv[1]) return("Clonal")
  return("Subclonal")
  }
