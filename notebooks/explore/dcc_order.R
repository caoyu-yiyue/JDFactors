suppressPackageStartupMessages({
  library(rmgarch)
  library(xts)
  library(foreach)
  library(doParallel)
  library(PerformanceAnalytics)
})

source("src/data/read_data.R")
source("src/config.R")
source("src/process/multi_garch_mdl.R")
source("src/process/copula_mdl.R")
source("src/process/rolling_copula.R")
source("src/process/weight_optimize.R")
source("src/result/port_ret.R")

facs_xts <- read_fac_xts()
in_sample_end <- in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS)
multigarchfit_list <- read_rolling_multigarchfit()
rolling_mean <- read_rolling_mean()

arma_order_for_roll <- matrix(rep(3, 10), nrow = 2)
multi_garch_spec <- all_facs_multigarch(
  arma_order_df = arma_order_for_roll,
  fit = FALSE
)

# c(1, 2)
t_dcc_12 <- fit_garch_copula(
  multigarch_spec = multi_garch_spec,
  copula_type = "mvt", is_dcc = TRUE,
  dcc_order = c(1, 2), fit = FALSE
)

t_dcc_12_rcov <- rolling_cgarch_rcov(
  data = facs_xts,
  pure_cgarch_spec = t_dcc_12,
  start_row = in_sample_end,
  step_by = ROLLING_STEP,
  multigarchfit_list = multigarchfit_list
)

# c(2, 1)
t_dcc_21 <- fit_garch_copula(
  multigarch_spec = multi_garch_spec,
  copula_type = "mvt", is_dcc = TRUE,
  dcc_order = c(2, 1), fit = FALSE
)

t_dcc_21_rcov <- rolling_cgarch_rcov(
  data = facs_xts,
  pure_cgarch_spec = t_dcc_21,
  start_row = in_sample_end,
  step_by = ROLLING_STEP,
  multigarchfit_list = multigarchfit_list
)

# c(2, 2)
t_dcc_22 <- fit_garch_copula(
  multigarch_spec = multi_garch_spec,
  copula_type = "mvt", is_dcc = TRUE,
  dcc_order = c(2, 2), fit = FALSE
)

t_dcc_22_rcov <- rolling_cgarch_rcov(
  data = facs_xts,
  pure_cgarch_spec = t_dcc_22,
  start_row = in_sample_end,
  step_by = ROLLING_STEP,
  multigarchfit_list = multigarchfit_list
)

t_dcc_list <- list(
  t_dcc_12 = t_dcc_12_rcov, t_dcc_21 = t_dcc_21_rcov,
  t_dcc_22 = t_dcc_22_rcov
)

# 解最优化
rcov_names <- names(t_dcc_list)
n_fac <- ncol(facs_xts)
gammas <- c(3, 8, 20, 50)
opt_weights <- list()
for (gamma in gammas) {
  for (name in rcov_names) {
    rcov_xts <- t_dcc_list[[name]]
    mean_rcov_merged <- cbind(rolling_mean, rcov_xts)
    opt_weight <- rolling_opt(
      data = mean_rcov_merged, gamma = gamma,
      n_fac = n_fac, sum_1 = TRUE
    )
    # 最优化过程中小于0 的特别小的值，直接变成0
    round_0 <- replace(opt_weight, opt_weight < 0, 0)
    opt_weights[[as.character(gamma)]][[name]] <- round_0
  }
}

# 求组合收益率
port_rets_list <- rets_for_all_sigma(
  opt_weights = opt_weights,
  facs_ret = facs_xts
)

basic_static_table <- function(ret_xts, data_freq) {
  scale <- switch(data_freq, "Week" = 52, "Day" = 252, "Month" = 12)
  annual_table <- table.AnnualizedReturns(ret_xts, scale = scale)
  return(annual_table)
}

basic_static_table(ret_xts = port_rets_list[["50"]], data_freq = "Week")
