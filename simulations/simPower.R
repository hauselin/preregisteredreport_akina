setwd("/scratch/m/mickeyi/hauselin/AkinaPower/")

library(data.table); library(paramtest); library(lme4); library(lmerTest); library(tidyverse); library(parallel)

mlm_test <- function(simNum, N, b1, b2, b0 = 0, varInt = 1, varSlope_b1 = 1, varSlope_b2 = 1, varResid = 1) {
  # timePoints <- 4

  b1_levels <- 4
  b2_levels <- 5
  rows_per_subject <- b1_levels * b2_levels

  combinations <- expand.grid(x1 = 0:(b1_levels-1), x2 = 0:(b2_levels-1))

  dt1 <- data.table::data.table(subject = rep(1:N, each = nrow(combinations)), x1 = rep(combinations$x1, N), x2 = rep(combinations$x2, N))

  dt1[, sub_int := rep(rnorm(N, 0, sqrt(varInt)), each = nrow(combinations))] # random intercept for each subject
  dt1[, sub_slope_b1 := rep(rnorm(N, 0, sqrt(varSlope_b1)), each = nrow(combinations))] # random slope b1 for each subject
  dt1[, sub_slope_b2 := rep(rnorm(N, 0, sqrt(varSlope_b2)), each = nrow(combinations))] # random slope b2 for each subject

  # generate data
  # y-intercept as a function of b0 plus random intercept;
  # slope as a function of b1 plus random slope
  dt1[, y := (b0 + sub_int) + (b1 + sub_slope_b1) * x1 + (b2 + sub_slope_b2) * x2 + rnorm(.N, 0, sqrt(varResid))]

  # for more complex models that might not converge, tryCatch() is probably
  # a good idea
  return <- tryCatch({
    model <- lmerTest::lmer(y ~ x1 + x2 + (1 | subject), data=dt1)
    # when using parallel processing, we must refer to functions from
    # packages directly, e.g., package::function()

    # summary(model)
    est_int <- coef(summary(model))['(Intercept)', 'Estimate']
    est_x1 <- coef(summary(model))['x1', 'Estimate']
    est_x2 <- coef(summary(model))['x2', 'Estimate']
    se_x1 <- coef(summary(model))['x1', 'Std. Error']
    se_x2 <- coef(summary(model))['x2', 'Std. Error']

    df_x1 <- coef(summary(model))['x1', 'df']
    df_x2 <- coef(summary(model))['x2', 'df']
    t_x1 <- coef(summary(model))['x1', 't value']
    t_x2 <- coef(summary(model))['x2', 't value']
    r_x1 <- sqrt(t_x1^2 / (t_x1^2 + df_x1))
    r_x2 <- sqrt(t_x2^2 / (t_x2^2 + df_x2))

    p_x1 <- coef(summary(model))['x1', 'Pr(>|t|)']
    p_x2 <- coef(summary(model))['x2', 'Pr(>|t|)']

    return(c(est_int = est_int,
             est_x1 = est_x1, est_x2 = est_x2,
             se_x1 = se_x1, se_x2 = se_x2,
             df_x1 = df_x1, df_x2 = df_x2,
             t_x1 = t_x1, t_x2 = t_x2,
             r_x1 = r_x1, r_x2 = r_x2,
             sig_x1 = p_x1 < 0.02,
             sig_x2 = p_x2 < 0.02))
  },
  error=function(e) {
    #message(e)  # print error message
    return(c(est_int = NA,
             est_x1 = NA, est_x2 = NA,
             se_x1 = NA, se_x2 = NA,
             df_x1 = NA, df_x2 = NA,
             t_x1 = NA, t_x2 = NA,
             r_x1 = NA, r_x2 = NA,
             sig_x1 = NA,
             sig_x2 = NA))
  })

  return(return)
}

power_mlm <- grid_search(mlm_test, 
params <- list(N = c(60, 70, 80, 90), b1 = c(0.05, 0.08, 0.12), b2 = c(0.05, 0.08, 0.12),
                                                 varInt = c(0.01, 0.05, 0.1), 
                                                 varSlope_b1 = c(0.01, 0.05, 0.1), varSlope_b2 = c(0.01, 0.05, 0.1), 
                                                 varResid = c(0.01, 0.1, 0.2, 0.3)),
                         n.iter = 1000, output = 'data.frame', 
                         parallel = "snow", ncpus = detectCores())

print('iterations')
print(nrow(power_mlm$results))
print(power_mlm$timing)

print('seconds per iteration')
x1 <- power_mlm$timing[3] / nrow(power_mlm$results)
print(x1)

print('hours for 5000000 iterations')
print(x1 * 5000000 / 60 / 60)

power_mlm_tbl <- data.table(power_mlm$results)
# print(power_mlm_tbl)

write_rds(power_mlm, "power_simulations.RDS")
fwrite(power_mlm_tbl, "power_simulations.csv")

