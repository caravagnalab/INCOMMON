##
x = readRDS("./testdata.rds")
test = TAPACLOTH::run_classifier(x, alpha_level = 0.05, rho = 0.01, model = "beta-binomial")
model = "beta-binomial"
gene_name = "SPEN"
alpha_level = 0.05
rho = 0.01

plot_test = function(x){
  
  plotmodels = lapply(names(x$classifier), function(model){
    
    plotlist = lapply(x$classifier[[model]]$data$gene %>% unique(), function(g){
      
      nvtest = x$classifier[[model]]$data %>% 
        filter(gene == g) %>% 
        pull(NV) %>% 
        unique()
      
      purity = x$purity
      
      gdata = x$classifier[[model]]$data %>% dplyr::filter(gene == g)
      
      y = lapply(1:(
        gdata %>% nrow()
      ),
      function(i) {
        
        dp = gdata[i,]$DP
        k = gdata[i,]$karyotype
        m = gdata[i,]$multiplicity
        ploidy = stringr::str_split(k, pattern = ":")[[1]] %>% as.integer() %>% sum()
        
        if ((model %>% tolower()) == "binomial") {
          p = dbinom(
            x = 1:dp,
            size = dp,
            prob = m * purity / (2 * (1 - purity) + purity * ploidy)
          )
        }
        if ((model %>% tolower()) == "beta-binomial") {
          p = VGAM::dbetabinom(
            x = 1:dp,
            size = dp,
            rho = rho,
            prob = m * purity / (2 * (1 - purity) + purity * ploidy)
          )
        }
        tibble(
          nv = 1:dp,
          p = p,
          class = gdata[i,]$class,
          outcome = ifelse(gdata[i,]$pvalue > alpha_level, "PASS", "FAIL"),
          l_a = gdata[i,]$l_a,
          r_a = gdata[i,]$r_a
        )
      }) %>% do.call(rbind, .)
      y = y %>% 
        mutate(ff = case_when(
          outcome == "FAIL" ~ ggplot2::alpha("gray", 0.4),
          outcome == "PASS" & class == "k=1:0,m=1" ~ "steelblue",
          outcome == "PASS" & class == "k=1:1,m=1" ~ ggplot2::alpha("forestgreen", 0.8),
          outcome == "PASS" & class == "k=2:0,m=1" ~ "turquoise4",
          outcome == "PASS" & class == "k=2:0,m=2" ~ ggplot2::alpha("turquoise4",0.5),
          outcome == "PASS" & class == "k=2:1,m=1" ~ ggplot2::alpha("orange",0.8),
          outcome == "PASS" & class == "k=2:1,m=2" ~ ggplot2::alpha("orange",0.5),
          outcome == "PASS" & class == "k=2:2,m=1" ~ "firebrick3",
          outcome == "PASS" & class == "k=2:2,m=2" ~ ggplot2::alpha("firebrick3",0.5)
        ))
      
      pp = lapply(arrange(y, outcome)$class %>% unique(), function(c){
        # ff = case_when(
        #   unique(filter(y, class==c)$outcome) == "FAIL" ~ ggplot2::alpha("gray", 0.4),
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=1:0,m=1" ~ "steelblue",
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=1:1,m=1" ~ ggplot2::alpha("forestgreen", 0.8),
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:0,m=1" ~ "turquoise4",
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:0,m=2" ~ ggplot2::alpha("turquoise4",0.5),
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:1,m=1" ~ ggplot2::alpha("orange",0.8),
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:1,m=2" ~ ggplot2::alpha("orange",0.5),
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:2,m=1" ~ "firebrick3",
        #   unique(filter(y, class==c)$outcome) == "PASS" & c == "k=2:2,m=2" ~ ggplot2::alpha("firebrick3",0.5)
        # )
        geom_area(
          y %>% filter(class==c,nv > l_a &
                         nv < r_a),
          mapping = aes(x = nv, y = p, fill = ff),
          position = 'identity'
        )
      })
      
      ggplot() + (pp %>% unlist()) +
        scale_fill_manual(
          name = "class",
          values = c(
            "#BEBEBE66" = "#BEBEBE66",
            "steelblue" = "steelblue",
            "#228B22CC" = "#228B22CC",
            "turquoise4" = "turquoise4",
            "#00868B80" = "#00868B80",
            "#FFA500CC" = "#FFA500CC",
            "#FFA50080" = "#FFA50080",
            "firebrick3" = "firebrick3",
            "#CD262680" = "#CD262680"
          ),
          labels = c(
            "FAIL",
            "k=1:0,m=1",
            "k=1:1,m=1",
            "k=2:0,m=1",
            "k=2:0,m=2",
            "k=2:1,m=1",
            "k=2:1,m=2",
            "k=2:2,m=1",
            "k=2:2,m=2"
          )
        ) + ggplot2::geom_vline(xintercept = nvtest,
                              linetype = "longdash",
                              color = "black")+
        CNAqc:::my_ggplot_theme() +
        ggplot2::labs(
          x = 'NV',
          y = "P(X = NV)",
          caption = paste0("Test using ", model, " model"),
          title = paste0("Gene ",g,": coverage ", max(gdata$DP[i]), ' with purity ', x$purity),
          subtitle = paste0('Alpha-level: ', alpha_level)
        )
    })
    names(plotlist) = x$classifier[[model]]$data$gene %>% unique()
    plotlist
  })
  names(plotmodels) = names(x$classifier)
  plotmodels
}

plot_test(x)


# ggplot2::ggplot() +
#   ggplot2::geom_vline(xintercept = nvtest,
#              linetype = "longdash",
#              color = "black") +
#   ggplot2::geom_line(y, mapping = aes(x = nv, y = p)) +
#   # ggplot2::geom_vline(y, mapping = aes(xintercept = l_a)) +
#   # ggplot2::geom_vline(y, mapping = aes(xintercept = r_a)) +
#   ggplot2::geom_area(y %>% filter(nv>l_a & nv<r_a), mapping = aes(x=nv, y=p), alpha=0.5)+
#   CNAqc:::my_ggplot_theme() +
#   ggplot2::facet_wrap( ~ class,ncol = 2)
#   # ggplot2::labs(
#   #   x = 'NV',
#   #   y = "P(X = NV)",
#   #   caption = model,
#   #   title = paste0("Coverage ", dp, ' with purity ', purity),
#   #   subtitle = paste0('Alpha-level: ', alpha_level))
# ##
# 
# ggplot2::ggplot() +
#   # ggplot2::geom_line(y, mapping = aes(x = nv, y = p, color = class), fill = NA) +
#   ggplot2::geom_area(
#     y %>% filter(nv > l_a &
#                    nv < r_a),
#     mapping = aes(x = nv, y = p, fill = class, alpha = outcome),
#     position = 'identity'
#   ) + scale_alpha_discrete(range=c(0.08,0.9))+
#   # scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, name = 'Set3'))+
#   scale_fill_manual(values = ggsci::pal_lancet(palette = "lanonc")(8))+
#   ggplot2::geom_vline(xintercept = nvtest,
#                       linetype = "longdash",
#                       color = "black") +
#   CNAqc:::my_ggplot_theme() 



