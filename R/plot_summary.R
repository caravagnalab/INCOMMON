#'Plotting function for class \code{'TAPACLOTH'}.
#' @description
#' Generates a composite heatmap with summary statistics of classification over multiple samples.
#' @param x A list of obj of class \code{'TAPACLOTH'} containing a `classifier`.
#' @param model Model used for the classification, either "Binomial" or "Beta-Binomial".
#' @import pheatmap
#' @import CNAqc
#' @import ggplot2
#' @import ggsci 
#' @import ggridges 
#' @importFrom dplyr filter mutate rename select %>% 
#' @return An object of class \code{'pheatmap'}.
#' @export
#' @examples
#' x = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity)
#' x = run_classifier(x)
#' plot_heatmap(x)
plot_summary = function(x,
                        model = 'binomial',
                        out_file = "./tapacloth_pheatmap.pdf",
                        main = "summary statistics of classification",
                        n_frequency_colors = 4
)
{
  # stopifnot(inherits(x, "TAPACLOTH"))
  
  model = model %>% tolower()
  
  # Tibble 2 matrix
  t2m = function(x){
    m =  x[, 2:ncol(x)] %>% as.matrix()
    colnames(m) = names(x)[2:length(names(x))]
    rownames(m) = x$gene
    m
  }
  
  ## Main matrix with classification results
  input_classes_matrix = x %>% 
    dplyr::filter(label!="out of sample") %>% 
    group_by(gene,label) %>% 
    summarise(n=n()) %>% 
    summarise(label,f=n/sum(n)) %>% 
    tidyr::pivot_wider(names_from = label, values_from = f, values_fill=0)
  
  classes_matrix = input_classes_matrix %>% t2m()
  
  ## Uncertainty
  input_uncertainty_matrix = x %>% 
    group_by(gene) %>% 
    summarise(uncertainty = median(uncertainty))
  
  uncertainty_matrix = input_uncertainty_matrix %>% t2m()
  
  ## Occurrence
  input_occurrence_matrix = x %>% 
    group_by(gene) %>% 
    summarise(N = n())
  
  occurrence_matrix = input_occurrence_matrix %>% t2m()
  
  ## Gene role
  input_gene_role_matrix = x %>% 
    group_by(gene) %>% 
    summarise(gene_role = unique(gene_role))
  
  gene_role_matrix = input_gene_role_matrix %>% t2m()
  
  # roles_all = x$gene_role %>% unique()
  # roles_all = strsplit(roles_all, ',') %>% unlist %>% unique
  # roles_all_others = roles_all[!(roles_all %in% c("TSG", "oncogene", NA))]
  # 
  # df_roles = x %>%
  #   dplyr::select(gene, gene_role) %>%
  #   dplyr::rename(role = gene_role) %>%
  #   dplyr::distinct() %>%
  #   rowwise() %>%
  #   mutate(
  #     TSG = grepl('TSG', role),
  #     oncogene = grepl('oncogene', role)
  #   ) %>%
  #   dplyr::select(-role) %>%
  #   as.data.frame()
  # 
  # annotation_rows = df_roles
  # rownames(annotation_rows) = annotation_rows$gene
  # annotation_rows = annotation_rows[,2:3]
  # annotation_rows = annotation_rows[order(row.names(annotation_rows)),]
  # 
  # annotation_rows = apply(annotation_rows, 2, function(x)
  #   ifelse(x, "YES", "NO")) %>% data.frame
  # 
  # annotation_rows = cbind(annotation_rows, uncertainty_matrix)


  # Order rows (first TSG then oncogene) and columns (increasing ploidy)
  # classes_matrix = classes_matrix[c(which(annotation_rows$TSG == "YES"),
  #                                   which(annotation_rows$oncogene == "YES", )),
  #                                 order(colnames(classes_matrix))]
  
  annotation_rows = cbind(gene_role_matrix, uncertainty_matrix, occurrence_matrix) %>% 
    as.data.frame()
  annotation_rows[,"uncertainty"] = annotation_rows[,"uncertainty"] %>% as.numeric()
  annotation_rows[,"N"] = annotation_rows[,"N"] %>% as.numeric()
  
  classes_matrix = classes_matrix[c(which(annotation_rows[,"gene_role"] == "TSG"),
                                    which(annotation_rows[,"gene_role"] == "oncogene")), 
                                  order(colnames(classes_matrix))]
  
  presence_annotation = 'darkgray'
  onco_color = 'steelblue'
  tsg_color = 'indianred'
  
  # annotation_colors = list(
  #   uncertainty = hcl.colors(20,"Grays", rev = T),
  #   TSG = c(`NO` = "white", `YES` = tsg_color),
  #   oncogene = c(`NO` = "white", `YES` = onco_color)
  # )
  
  annotation_colors = list(
    uncertainty = hcl.colors(20,"Grays", rev = T),
    gene_role = c("TSG" = tsg_color, "oncogene" = onco_color),
    N = hcl.colors(20,"Greens 3", rev = T)
  )
  
  pheatmap::pheatmap(
    mat = classes_matrix,
    color = hcl.colors(n_frequency_colors,"RdPu", rev = T),
    main = main,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    fontsize_row = 6,
    fontsize = 7,
    cellwidth = 20,
    cellheight = 10,
    annotation_row = annotation_rows,
    annotation_colors = annotation_colors,
    filename = out_file
  ) 
  
}
