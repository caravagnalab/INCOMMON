plot_test = function(x, id, model){
  stopifnot(inherits(x,"TAPACLOTH"))
  colors = CNAqc:::get_karyotypes_colors(c('1:0', '1:1', '2:0', '2:1', '2:2'))
  names_colors = expand.grid(names(colors), 1:2) %>% apply(1, paste, collapse = ' ')
  colors = sapply(names_colors, function(n)
    colors[strsplit(n, ' ')[[1]][1]])
  names(colors) = names_colors
  
  colors = c(colors, `out_of_sample` = 'gray')
  
  mdata = x %>% get_classifier(model) %>%
    get_data() %>%
    filter(id == !!id)
  
  dataset = pull(mdata, density)[[1]]
  
  dataset %>%
    ggplot() +
    geom_line(aes(x = NV, y = density), color = 'gainsboro',  size = .3) +
    CNAqc:::my_ggplot_theme() +
    scale_color_manual(values = colors) +
    geom_point(data = dataset %>% maximise(),
               aes(x = NV,
                   y = density,
                   color = label),
               size = 1) +
    geom_line(data = dataset %>% maximise(),
              aes(x = NV,
                  y = density,
                  color = label),
              size = .3) +
    labs(
      title = paste0(
        mdata$ref,
        ">",
        mdata$alt,
        "; ",
        mdata$gene,
        " (",
        mdata$gene_role,
        ")"
      ),
      subtitle = paste0(
      "Alleles: ",
      mdata$ploidy,
      "N (Mutated: ",
      mdata$multiplicity,
      "N, Wild Type: ",
      mdata$wt,
      "N)\nPurity: ",
      get_purity(x),
      ' - cutoff: ',
      x %>% get_params(model = model) %>% pull(cutoff)
    ),
    caption = paste0("Classified using ", model, " model")) +
    # guides(color = 'none') +
    # geom_hline(data = dataset %>% group_by(label) %>% summarise(cutoff) %>% unique() %>% filter(label != "out of sample"),
    #     mapping = aes(yintercept = cutoff,
    #     color = label),
    #     linetype = 'dashed',
    #     size = .5
    #   )
    geom_vline(
      xintercept = get_NV(x, id),
      color = ifelse(is.na(mdata$ploidy), 'indianred3', 'forestgreen'),
      linetype = 'dashed',
      size = .5
    )
}