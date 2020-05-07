suppressPackageStartupMessages({
  library(rmgarch)
  library(xts)
  library(foreach)
})

source("src/data/read_data.R")
source("src/process/multi_garch_mdl.R")
source("src/process/copula_mdl.R")


set_cgarchspec_fixed <- function(cspec_obj, cfit_obj) {
  #' @title 为cgarchspec 对象设定fixed parameter
  #'
  #' @param cspec_obj cGARCHspec 对象。一个等待fix param 的spec
  #' @param cfit_obj cGARCHfit 对象。一个用于获取param 的对象
  #' @details 本质上是将cfit_obj 中的参数赋值到cspec_obj 中。
  #' 未来可能会允许外部传入自定义参数。
  #' 算法来自https://github.com/cran/rmgarch/blob/master/inst/
  #' rmgarch.tests/rmgarch.test3.R, Line 404
  #' @return cGARCHspec 对象。设定完fixed param 的spec。

  # 首先对每个单变量GARCH 模型fix parameter
  fac_num <- length(cspec_obj@umodel$fixed.pars)
  for (i in 1:fac_num) {
    cspec_obj@umodel$fixed.pars[[i]] <-
      as.list(cfit_obj@model$mpars[cfit_obj@model$midx[, i] == 1, i])
  }
  # 对cGARCHspec 对象固定Joint parameter
  setfixed(cspec_obj) <-
    as.list(cfit_obj@model$mpars[
      cfit_obj@model$midx[, (fac_num + 1)] == 1, (fac_num + 1)
    ])
  return(cspec_obj)
}


# 1. data 部分
facs_xts <- read_fac_xts()
year_endponits <- endpoints(facs_xts, on = "year")[-1]
sixth_year_endpoint <- year_endponits[6]


# 2. 指定cGARCHspec 部分
arma_order_for_roll <- matrix(rep(3, 10), nrow = 2)
multi_garch_spec <- all_facs_multigarch(
  arma_order_df = arma_order_for_roll,
  fit = FALSE
)
cop_garch_spec <- fit_garch_copula(
  multigarch_spec = multi_garch_spec,
  copula_type = "mvt", is_dcc = TRUE, fit = FALSE
)

# ================= rolling funciton part ================= #
full_rows <- nrow(facs_xts)
fit_time <- seq.int(
  from = sixth_year_endpoint,
  to = full_rows, by = 12
)
cls <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cls, cores = 2)
rolling_rcov <- foreach(t = fit_time, .combine = "rbind") %dopar% {
  # 3. fit 一个garch-copula 模型
  current_fit <- cgarchfit(spec = cop_garch_spec, data = facs_xts[1:t, ])

  # 4. 设定fixed param 并filter 部分
  cspec_fixed_param <- set_cgarchspec_fixed(
    cspec_obj = cop_garch_spec,
    cfit_obj = current_fit
  )

  # 当总行数与t 相比大于等于12，则filter data 为到t 往后12期，预测期为12
  # 否则，filter data 为整个facs_xts，预测期为所有数据行 - t
  if (full_rows - t >= 12) {
    filter_data <- facs_xts[1:(t + 12), ]
    forcast_t <- 12
  } else {
    filter_data <- facs_xts
    forcast_t <- full_rows - t
  }

  current_filter <- cgarchfilter(cspec_fixed_param,
    data = filter_data, # 这里在数据的最末尾会出现out of bound 问题
    filter.control = list(n.old = t)
  )

  # 5. 返回需要的部分
  forcasted_cov <- tail(rcov(current_filter, output = "matrix"), forcast_t)
  return(forcasted_cov)
}
stopCluster(cls)
