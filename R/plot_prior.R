#' Visualize prior distribution for a gene (tumor-specific or pancancer).
#'
#' @param x A prior distribution in the format required for \code{INCOMMON},
#' such as \code{INCOMMON::pcawg_priors}.
#' @param tumor_type Tumor type for tumor-specific prior ('PANCA' for pan-cancer).
#' @param gene Gene for gene-specific prior.
#' @return An object or a list of objects of class \code{'ggplot2'}.
#' @export
#' @examples
#' # First load example classified data
#' data(MSK_classified)
#' # Plot classification results for a specific sample
#' plot_prior(x = MSK_classified, gene = 'TP53', tumor_type = 'PAAD')
#' @importFrom dplyr filter mutate rename select %>%
plot_prior = function(x, gene, tumor_type){

  get_ploidy <- function(label) {
    numbers <- unlist(strsplit(label, " and "))
    numbers <- as.numeric(numbers)
    total <- sum(numbers, na.rm = TRUE)
    return(total)
  }

  get_prior(x, gene, tumor_type) %>%
    dplyr::mutate(
      state = dplyr::case_when(
        label %in% c("2N (Mutated: 1N)") ~ "HMD",
        label %in% c("4N (Mutated: 1N)", "3N (Mutated: 1N)") ~ "Tier-2",
        label %in% c("1N (Mutated: 1N)") ~ "LOH",
        label %in% c("2N (Mutated: 2N)") ~ "CNLOH",
        label %in% c("3N (Mutated: 2N)", "4N (Mutated: 2N)") ~ "AM"
      )
    ) %>%
    ggplot2::ggplot()+
    ggplot2::geom_bar(ggplot2::aes(
      x = '',
      y = p,
      fill = state
    ), stat = 'identity')+
    scale_color_INCOMMON_class()+
    ggplot2::coord_flip()+
    CNAqc:::my_ggplot_theme(cex = .8)+
    ggplot2::xlab('')+ggplot2::ylab('')+
    ggplot2::facet_wrap(~tumor_type~gene)+
    ggplot2::guides(fill = ggplot2::guide_legend(title = 'INCOMMON class'), ncol = 2)
}
