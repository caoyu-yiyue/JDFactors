#!/usr/bin/env Rscript

#############################################################################
# 计算滚动的multigarch fit 对象，成为下一步rolling_cop 的基础
#############################################################################

suppressPackageStartupMessages({
  library(rugarch)
  library(foreach)
  library(parallel)
  library(doParallel)
})

source("src/data/read_data.R")
source("src/config.R")
source("src/process/multi_garch_mdl.R")


rolling_multigarch_fit <- function(data, multigarch_spec, start_t, step_by) {
  #' @title 滚动计算multigarch fit 对象，然后返回这些对象的list
  #'
  #' @param data xts 对象。进行rolling fit 的因子数据
  #' @param multigarch_spec rugarch::multigarchSpec 对象。滚动时共用的spec
  #' @param start_t int. 最开始估计时的起始数据行数
  #' @param step_by int. 每隔多久refit 一次
  #' @return list of uGARCHmultifit. 滚动估计的fit 对象们，name 为所用数据的行数

  # 根据传入的开始时间和step，计算出需要refit 的行
  fit_time <- seq.int(
    from = start_t,
    to = nrow(data),
    by = step_by
  )

  # 并行计算multi garch fit，所有的fit 对象组合返回一个list
  cls <- parallel::makeForkCluster(parallel::detectCores())
  doParallel::registerDoParallel(cls)
  rolling_multigarch_fits <- foreach(t = fit_time) %dopar% {
    multi_garch_fit <- multifit(
      multigarch_spec,
      data = tryCatch(data[1:t, ], error = function() data),
      solver = "solnp"
    )
    return(multi_garch_fit)
  }
  parallel::stopCluster(cls)

  # 对list 命名并返回
  names(rolling_multigarch_fits) <- fit_time
  return(rolling_multigarch_fits)
}


rolling_multigarch_main <- function() {
  #' @title rolling 计算multigarch fit 的主函数，以此为基础来计算后面不同的copula 参数

  # 读取数据，并找到起始行
  facs_xts <- read_fac_xts()
  in_sample_end <- in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS)

  # 设定每次refit 共用的multigarch spec 对象
  arma_order_for_roll <- matrix(rep(3, 10), nrow = 2)
  multigarch_spec <- all_facs_multigarch(arma_order_for_roll, fit = FALSE)

  rolling_multigarch_fits <- rolling_multigarch_fit(
    data = facs_xts, multigarch_spec = multigarch_spec,
    start_t = in_sample_end, step_by = ROLLING_STEP
  )

  # 为list 添加名字，每个项目名字为该fit_time(实际是行数)。最后保存。
  args <- commandArgs(trailingOnly = TRUE)
  saveRDS(object = rolling_multigarch_fits, file = args[[1]])
}


if (!interactive()) {
  rolling_multigarch_main()
}