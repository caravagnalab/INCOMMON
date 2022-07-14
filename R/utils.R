get_data = function(x){
  return(x$data)
}

get_sample = function(x){
  return(x$sample)
}

get_purity = function(x){
  return(x$purity)
}

idify = function(x){
  y = get_data(x)
  y = y %>% 
    mutate(id = paste(chr,from,to,ref,alt,sep = ":"))
  x$data = y
  return(x)
}

unidify = function(x){
  y = get_data(x)
  y = y %>% 
    dplyr::select(-id)
  x$data = y
  return(x)
}

gene_data = function(x, gene_id){
  y = x$data %>% 
    dplyr::filter(gene == gene_id)
  return(y)
}

get_DP = function(x, mutation_id){
  x = idify(x)
  get_data(x) %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(DP)
}

get_NV = function(x, mutation_id){
  x = idify(x)
  get_data(x) %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(NV)
}

get_VAF = function(x, mutation_id){
  x = idify(x)
  get_data(x) %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(VAF)
}

get_gene = function(x, mutation_id){
  x = idify(x)
  x$data %>% 
    dplyr::filter(id == mutation_id) %>% 
    pull(gene)
}

get_pvalues = function(x, null_model, mutation_id){
  y = null_model$test %>% 
    dplyr::select(karyotype, multiplicity, l_a, r_a)
  
  y$pvalue = sapply(null_model$test$inputs, function(s) {
    s$p[get_NV(x, mutation_id)]
  })
  
  y$gene = get_gene(x, mutation_id)
  
  return(y)
}

plot_gene = function(x,model,gene_name){
  model = model %>% tolower()
  x$classifier[[model]]$plot_test[[gene_name]]
}

get_purity_estimate = function(x, model){
  y = x$purity_estimate[[model]]
  return(y)
}
get_purity_bmix = function(x, model){
  y = get_purity_estimate(x, model)
  return(y$purity)
}

get_reliability = function(x, model){
  y = get_purity_estimate(x, model)
  return(y$reliability)
}