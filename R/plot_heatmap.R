#'Plotting function for class \code{'TAPACLOTH'}.
#' @description
#' Generates a composite heatmap with useful information on data and classification results.
#' @param x An obj of class \code{'TAPACLOTH'} containing a `classifier`.
#' @import pheatmap
#' @import CNAqc
#' @import ggplot2
#' @import ggsci 
#' @import ggridges 
#' @importFrom dplyr filter mutate rename select %>% 
#' @return An object of class \code{'TAPACLOTH'} containing a list of ggplot2
#' plots named `plot_test` inside `classifier`.
#' @export
#' @examples
#' x = init(mutations = example_data$data, sample = example_data$sample, purity = example_data$purity)
#' x = run_classifier(x)
#' plot_heatmap(x)
plot_heatmap = function(x, model = 'binomial')
{
  stopifnot(inherits(x, "TAPACLOTH"))
  
  # Tibble 2 matrix
  t2m = function(x){
    m =  x[, 2:ncol(x)] %>% as.matrix()
    colnames(m) = names(x)[2:length(names(x))]
    rownames(m) = x$label
    m
  }
  
  # cli::cli_alert("Plotting model pheatmap")
  
  input_classes = x %>%
    TAPACLOTH:::get_classifier(model = model) %>%
    TAPACLOTH:::get_data() %>%
    dplyr::mutate(
      label = paste0(gene, " (", from, ") ", ref, '>', alt),
      label_k = paste(karyotype, ' (', multiplicity, ')', sep = '')
    )
  
  # Matrix of the classifications
  input_classes_matrix = input_classes %>%
    dplyr::mutate(outcome = ifelse(outcome, 1, 0)) %>%
    dplyr::select(label, label_k, outcome) %>%
    tidyr::pivot_wider(names_from = label_k, values_from = outcome)
  
  # Removed
  # input_classes_matrix$other = 0
  # no_classes = apply(input_classes_matrix[2:ncol(input_classes_matrix)], 1, function(x) all(x == 0))
  # input_classes_matrix$other[no_classes] = 1
  
  # Idiotic formats
  matrix_outcomes = input_classes_matrix %>% t2m()
  
  # p-values
  pvalue_matrix = input_classes %>%
    dplyr::select(label, label_k, pvalue) %>%
    tidyr::pivot_wider(names_from = label_k, values_from = pvalue)
  
  p_matrix =  pvalue_matrix %>% t2m
  
  # Nicola - I don't understand why I have to use <, shouldn't it be >= to remove non-significant ones?
  alpha_test = TAPACLOTH::get_params(x) %>% 
    dplyr::filter(model == !!model) %>% 
    pull(alpha)
  
  p_matrix[p_matrix < alpha_test] = NA
  
  p_matrix = apply(p_matrix, c(1, 2), function(x) {
    if (is.na(x)) return('')
    if (x < 0.33) return("*")
    if (x > 0.33 & x < 0.67) return("**")
    return("***")
  })
  
  # Read counts data
  df_counts = input_classes %>%
    dplyr::select(label, NV, DP, VAF) %>%
    dplyr::distinct() %>%
    data.frame()
  
  rownames(df_counts) = df_counts$label
  df_counts = df_counts[, 2:ncol(df_counts)]
  
  # Roles
  roles_all = input_classes$gene_role %>% unique()
  roles_all = strsplit(roles_all, ',') %>% unlist %>% unique
  roles_all_others = roles_all[!(roles_all %in% c("TSG", "oncogene", NA))]
  
  df_roles = input_classes %>%
    dplyr::select(label, gene_role) %>%
    dplyr::rename(role = gene_role) %>%
    dplyr::distinct() %>%
    rowwise() %>%
    mutate(
      TSG = grepl('TSG', role),
      oncogene = grepl('oncogene', role),
      other = sapply(roles_all_others, function(x)
        grepl(x, role)) %>% any,
      unknown = is.na(role)
    ) %>%
    dplyr::select(-role) %>%
    as.data.frame()
  
  rownames(df_roles) = df_roles$label
  df_roles = df_roles[,-1]
  
  # df_roles = apply(df_roles, 2, function(x) ifelse(x, 1, 0)) %>% data.frame
  df_roles = apply(df_roles, 2, function(x)
    ifelse(x, "YES", "NO")) %>% data.frame
  
  # Annotation rows
  annotation_rows = cbind(df_counts, df_roles)
  
  # I decided we do not show these
  annotation_rows = annotation_rows[,!apply(annotation_rows, 2, function(x)
    all(x == "NO"))]
  annotation_rows = annotation_rows[, colnames(annotation_rows) != 'NV']
  annotation_rows = annotation_rows[, colnames(annotation_rows) != 'other']
  
  presence_annotation = 'darkgray'
  onco_color = 'orange'
  tsg_color = 'steelblue'
  
  annotation_colors = list(
    TSG = c(`NO` = "white", `YES` = tsg_color),
    oncogene = c(`NO` = "white", `YES` = onco_color),
    other = c(`NO` = "white", `YES` = presence_annotation),
    unknown = c(`NO` = "white", `YES` = presence_annotation),
    # DP = RColorBrewer::brewer.pal(10, 'Purples'),
    DP = viridis::viridis(10, direction = -1),
    NV = RColorBrewer::brewer.pal(8, 'Purples'),
    # VAF = RColorBrewer::brewer.pal(10, 'Spectral')
    VAF = viridis::viridis(10, direction = -1)
  )
  
  pheatmap::pheatmap(
    mat = matrix_outcomes,
    main = x$sample,
    cluster_cols = FALSE,
    fontsize_row = 6,
    fontsize = 7,
    cellwidth = 20,
    cellheight = 10,
    display_numbers = p_matrix,
    number_color = 'white',
    # border_color = NA
    gaps_col = c(1, 2, 4, 6, 8),
    breaks = c(-1, 0, 1),
    color = c(`0` = ggplot2::alpha('gainsboro', .5), `1` = 'indianred3'),
    legend_breaks = c(0, 1),
    legend_labels = c(`0` = "Reject H0", `1` = "Fail to\nReject H0"),
    annotation_row = annotation_rows,
    annotation_colors = annotation_colors
  )
  
}
