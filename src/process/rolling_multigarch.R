#!/usr/bin/env Rscript

#############################################################################
# 计算滚动的multigarch fit 对象，成为下一步rolling_cop 的基础
#############################################################################

suppressPackageStartupMessages({
  library(rugarch)
  library(foreach)
  library(parallel)
  library(doParallel)
  library(optparse)
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
    # 首先使用multifit，建立一个multiGARCHfit obj
    multi_garch_fit <- multifit(
      multigarch_spec,
      data = tryCatch(data[1:t, ], error = function(cond) data),
      # fit.control = list(scale = 10000),
      solver = "hybrid"
    )

    # 判定是否每个变量都解开了最优化，是否有cvar，如果没有，则对该变量单独fit
    conver_flags <- .conver_for_multigarchfit(multi_garch_fit)
    cvar_flags <- .cvar_for_multigarchfit(multi_garch_fit)
    problem_idx <- which(!(conver_flags & cvar_flags))
    for (i in problem_idx) {
      # 每个fac 尝试最多5 次
      for (try_time in 1:5) {
        ugfit <- ugarchfit(
          spec = multigarch_spec@spec[[i]],
          data = tryCatch(data[1:t, i], error = function(cond) data[, i]),
          solver = "hybrid", fit.control = list(scale = 10**try_time)
        )

        if (ugfit@fit$convergence == 0 & "cvar" %in% names(ugfit@fit)) {
          # 通过验证，被ugfit 放进multi garch fit 中
          multi_garch_fit@fit[[i]] <- ugfit
          break
        }

        # 如果到了最后一次还没有成功，尝试solver：lbfgs 五次
        if (try_time == 5) {
          for (lbfgs_try_time in 1:5) {
            ugfit <- ugarchfit(
              spec = multigarch_spec@spec[[i]],
              data = tryCatch(data[1:t, i], error = function(cond) data[, i]),
              solver = "lbfgs", fit.control = list(scale = 10**lbfgs_try_time)
            )
            if (ugfit@fit$convergence == 0 & "cvar" %in% names(ugfit@fit)) {
              # 通过验证，被ugfit 放进multi garch fit 中
              multi_garch_fit@fit[[i]] <- ugfit
              break
            }
          }
        }
      }
    }

    return(multi_garch_fit)
  }
  parallel::stopCluster(cls)

  # 对list 命名并返回
  names(rolling_multigarch_fits) <- fit_time
  return(rolling_multigarch_fits)
}


rolling_multigarch_main <- function() {
  #' @title rolling 计算multigarch fit 的主函数，以此为基础来计算后面不同的copula 参数

  # 解析命令行参数
  option_list <- list(
    make_option(
      opt_str = c("-f", "--data_freq"), type = "character",
      default = "Week", help = "Which data freq of factors?",
      metavar = "character"
    ),
    make_option(c("-o", "--out_put"),
      type = "character",
      default = NULL, help = "Path to save out put file.",
      metavar = "character"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opts <- optparse::parse_args(opt_parser)
  data_freq <- opts[["data_freq"]]

  # 读取数据，并找到起始行
  facs_xts <- read_fac_xts(data_freq = data_freq)
  in_sample_end <- in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS)

  # 设定每次refit 共用的multigarch spec 对象
  arma_order_for_roll <- matrix(rep(3, 10), nrow = 2)
  multigarch_spec <- all_facs_multigarch(arma_order_for_roll, fit = FALSE)

  rolling_multigarch_fits <- rolling_multigarch_fit(
    data = facs_xts, multigarch_spec = multigarch_spec,
    start_t = in_sample_end, step_by = ROLLING_STEP[data_freq]
  )

  # 验证rolling garch 全部fit成功
  source("test/rolling_garch_cvar_test.R")
  conv_result <- sapply(rolling_multigarch_fits, test_conv_cvar)
  if (FALSE %in% conv_result) {
    problem_idx <- which(!conv_result)
    print(problem_idx)
    stop("Not conv problem.")
  }

  # 保存生成的对象。
  saveRDS(object = rolling_multigarch_fits, file = opts[["out_put"]])
}


if (sys.nframe() == 0) {
  rolling_multigarch_main()
}
