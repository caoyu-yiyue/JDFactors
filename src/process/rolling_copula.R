suppressPackageStartupMessages({
  library(rmgarch)
  library(xts)
  library(foreach)
  library(doParallel)
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
rolling_cgarch_rcov <- function(data, pure_cgarch_spec,
                                start_row = 296, step_by = 12) {
  #' @title 对data 根据pure_cgarch_spec, 从第start_row 开始，每step_by(12期) refit
  #' 一次，保持参数固定向前filter 12 期，然后重复rolling fit，最终输出所有的预测
  #' cov 矩阵。
  #'
  #' @param data xts 对象，即多因子的xts 对象
  #' @param pure_cgarch_spec rmgarch::cgarchSpec。初始的cgarchSpec 对象。在这里
  #' 设定norm / t copula；static / dcc copula；以及指定固定参数
  #' @param start_row int, default 296. 最初进行cgarchfit 所用的数据行数，即in sample data。
  #' 默认是facs_xts 第6 年末的行。
  #' @param step_by int, default 12. 每隔多久重新fit 模型，默认12 期。
  #' @return xts 对象。返回的是每次fit 后向前filter step_by 期的方差-协方差矩阵rcov

  # 数据的总行数，以及需要refit 的行index 数。
  total_rows <- nrow(data)
  fit_time <- seq.int(
    from = start_row,
    to = total_rows, by = step_by
  )

  cls <- parallel::makeForkCluster(detectCores())
  doParallel::registerDoParallel(cls)
  rolling_rcov <- foreach(t = fit_time, .combine = "rbind") %dopar% {
    # 1. fit 一个garch-copula 模型
    current_fit <- cgarchfit(spec = pure_cgarch_spec, data = data[1:t, ])

    # 2. 设定fixed param 并filter 部分
    cspec_fixed_param <- set_cgarchspec_fixed(
      cspec_obj = pure_cgarch_spec,
      cfit_obj = current_fit
    )

    # 当总行数与t 相比大于等于step_by(12)，则filter data 为到t 往后step_by(12)期
    # 预测期为step_by(12)期；否则，filter data 为整个facs_xts，预测期为所有数据行 - t
    if (total_rows - t >= step_by) {
      filter_data <- data[1:(t + step_by), ]
      forcast_t <- step_by
    } else {
      filter_data <- data
      forcast_t <- total_rows - t
    }

    # 3. 根据固定参数的garch spec 向前filter
    current_filter <- cgarchfilter(cspec_fixed_param,
      data = filter_data, # 这里在数据的最末尾会出现out of bound 问题
      filter.control = list(n.old = t)
    )

    # 4. 返回需要的尾部部分
    forcasted_cov <- tail(rcov(current_filter, output = "matrix"), forcast_t)
    return(forcasted_cov)
  }
  parallel::stopCluster(cls)

  return(rolling_rcov)
}

# test_rcov <- rolling_cgarch_rcov(
#   data = facs_xts[1:(sixth_year_endpoint + 12 * 4)],
#   pure_cgarch_spec = cop_garch_spec
# )
