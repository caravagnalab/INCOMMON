plot_test_power = function(null_model)

{
  inputs = null_model$density
  coverage = null_model$coverage
  p = null_model$purity/2
  alpha_level = null_model$alpha_level
  nvs = null_model$density$nv %>% as.vector()
  l_a = null_model$nv[1]
  r_a = null_model$nv[2]
  l_v = null_model$vaf[1]
  r_v = null_model$vaf[2]

  model_string = "Binomial"
  if(null_model$model == 'Beta-Binomial')
    model_string = paste0("Beta-Binomial (rho = ", null_model$rho, ')')

  col_loh = 'steelblue'
  col_subclonal = 'indianred3'
  col_clonal = 'forestgreen'

  ggplot() +
    CNAqc:::my_ggplot_theme() +
    geom_rect(
      data = data.frame(
        xmin = 0,
        xmax = nvs[l_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_subclonal,
      alpha = .2
    ) +
    geom_rect(
      data = data.frame(
        xmax = coverage,
        xmin = nvs[r_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_loh,
      alpha = .2
    ) +
    geom_rect(
      data = data.frame(
        xmax = nvs[l_a],
        xmin = nvs[r_a],
        ymin = 0,
        ymax = 1
      ),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      fill = col_clonal,
      alpha = .2
    ) +
    geom_point(
      data = inputs %>% filter(nv < nvs[l_a]),
      aes(x = nv, y = p),
      size = .6,
      color = col_subclonal
    ) +
    geom_point(
      data = inputs %>% filter(nv > nvs[r_a]),
      aes(x = nv, y = p),
      size = .6,
      color = col_loh
    ) +
    geom_point(
      data = inputs %>% filter(nv >= nvs[l_a], nv <= nvs[r_a]),
      aes(x = nv, y = p),
      size = 1,
      color = col_clonal
    ) +
    geom_point(
      data = inputs,
      aes(x = nv, y = VAF),
      size = .6,
      shape = 3,
      color = 'gray'
    ) +
    geom_point(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(x = nv, y = VAF),
      size = 3,
      color = col_subclonal
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(
        x = nvs[l_a],
        y = VAF,
        xend = coverage,
        yend = VAF
      ),
      linetype = 'dashed',
      color = col_subclonal
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[l_a]),
      aes(
        x = nvs[l_a],
        y = VAF,
        xend = nvs[l_a],
        yend = 0
      ),
      linetype = 'dashed',
      color = col_subclonal
    ) +
    geom_point(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(x = nv, y = VAF),
      size = 3,
      color = col_loh
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(
        x = nvs[r_a],
        y = VAF,
        xend = coverage,
        yend = VAF
      ),
      linetype = 'dashed',
      color = col_loh
    ) +
    geom_segment(
      data = inputs %>% filter(nv == nvs[r_a]),
      aes(
        x = nvs[r_a],
        y = VAF,
        xend = nvs[r_a],
        yend = 1
      ),
      linetype = 'dashed',
      color = col_loh
    ) +
    geom_text(
      x = nvs[l_a] - 5,
      y = 1,
      label = 'subclonal',
      hjust = 1,
      size = 3,
      color = col_subclonal
    ) +
    geom_text(
      x = nvs[r_a] + 5,
      y = 0,
      label = 'clonal LOH',
      hjust = 0,
      size = 3,
      color = col_loh
    ) +
    labs(
      x = 'NV',
      y = "1 - P(X > NV)",
      caption = model_string,
      title = paste0("Coverage ", coverage, ' with purity ', p*2),
      subtitle = paste0('Alpha-level: ', alpha_level)
    ) +
    coord_cartesian(clip = 'off') +
    scale_y_continuous(sec.axis = sec_axis(~ ., name = 'VAF | coverage')) +
    theme(axis.line.y.right = element_line(color = "gray")) +
    annotation_custom(
      grid::textGrob(paste('>', round(r_v, 2)), gp = grid::gpar(col = col_loh, fontsize = 8)),
      xmin = coverage,
      xmax = coverage,
      ymin = r_v + .03,
      ymax = r_v + .03
    ) +
    annotation_custom(
      grid::textGrob(
        paste('<', round(l_v, 2)),
        gp = grid::gpar(col = col_subclonal, fontsize = 8)
      ),
      xmin = coverage,
      xmax = coverage,
      ymin = l_v - .03,
      ymax = l_v - .03
    ) +
    annotate(
      "text",
      x = 0,
      y = .95,
      label = "Subclonal",
      hjust = 0,
      size = 3,
      color = col_subclonal
    ) +
    annotate(
      "text",
      x = coverage,
      y = .95,
      label = "Clonal LOH",
      hjust = 1,
      size = 3,
      color = col_loh
    ) +
    annotate(
      "text",
      x = coverage * p,
      y = .95,
      label = "Clonal",
      hjust = .5,
      size = 3,
      color = col_clonal
    )
}

example_cartoons = function()
{
  ggpubr::ggarrange(
    plot_test_power(140, .3),
    plot_test_power(140, .3, model = 'betabinomial', rho = 0.001),
    plot_test_power(140, .3, model = 'betabinomial', rho = 0.01),
    ncol = 1,
    nrow = 3
  )

  ggpubr::ggarrange(
    plot_test_power(700, .4),
    plot_test_power(700, .4, model = 'betabinomial', rho = 0.001),
    plot_test_power(700, .4, model = 'betabinomial', rho = 0.01),
    ncol = 1,
    nrow = 3
  )
}
