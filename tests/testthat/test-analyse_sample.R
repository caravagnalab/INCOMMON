x = dplyr::tibble(sample = "test",
           gene = "test gene",
           nv = 50,
           dp = 100,
           vaf = 0.5,
           purity = 1
)

test_that("analyse_sample() returns a list", {
  expect_equal(analyse_sample(data = x,
                              sample_name = "test",
                              alpha_level = 1e-3,
                              model = "BetaBinomial",
                              rho = 0.01) %>% is.list(), TRUE)
})
