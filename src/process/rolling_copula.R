#!/usr/bin/env Rscript

#############################################################################
# 滚动计算copula 模型对应的rcov。
# 执行脚本时，将计算出t, norm 分别的dcc 和static copula 对应的rcov
# 需要传入一个参数：生成文件保存的位置
#############################################################################

suppressPackageStartupMessages({
  library(rmgarch)
  library(xts)
  library(foreach)
  library(doParallel)
})

source("src/data/read_data.R")
source("src/config.R")
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


multisol_cgarchfit <- function(spec, data, fit = NULL, ...) {
  #' @title 依此使用每个sovler 进行cgarchfit，顺序为"solnp", "nlminb", "gosolnp", "lbfgs"
  #'
  #' @param spec cGARCHspec 对象。
  #' @param data xts 对象。
  #' @param fit uGARCHmultifit 对象。已经拟合过的GARCH 模型
  #' @param ... 其他的参数，将传递给rmgarch::cgarchfit
  #' @details 函数将会依此使用"solnp", "nlminb", "gosolnp", "lbfgs" solver 求解，
  #' 直到返回一个cGARCHfit 对象。
  #' @return 如果拟合成功，将会返回一个rmgarch::cGARCHfit 对象
  #' 否则，如果拟合失败导致退化会返回一个rugarch::multiGARCHfit；
  #' 如果第四个solver 拟合失败，将会导致返回一个NULL 对象。

  # 依此使用每个solver 进行cgarch 拟合
  solvers <- c("solnp", "nlminb", "gosolnp", "lbfgs")
  for (solver in solvers) {
    # solver 使用时可能会出现错误，如果出错则直接赋值为NULL
    cgfit <- tryCatch(
      expr = {
        cgarchfit(
          spec = spec, data = data,
          fit = fit, solver = c("hybrid", solver), ...
        )
      },
      error = function(cond) {
        NULL
      }
    )

    # 如果class 是cGARCHfit 时，表示拟合成功，退出循环；否则会进入下一个循环
    if (is(cgfit, "cGARCHfit")) {
      break
    }

    # 当第四个solver 也失败时会运行到这里，抛出一个Warning。
    if (solver == solvers[length(solvers)]) {
      warning("Not convergence with every solver.")
    }
  }

  return(cgfit)
}


rolling_cgarch_rcov <- function(data, pure_cgarch_spec,
                                start_row, step_by,
                                multigarchfit_list = NULL) {
  #' @title 对data 根据pure_cgarch_spec, 从第start_row 开始，每step_by(12期) refit
  #' 一次，保持参数固定向前filter 12 期，然后重复rolling fit，最终输出所有的预测
  #' cov 矩阵。
  #'
  #' @param data xts 对象，即多因子的xts 对象
  #' @param pure_cgarch_spec rmgarch::cgarchSpec。初始的cgarchSpec 对象。在这里
  #' 设定norm / t copula；static / dcc copula；以及指定固定参数
  #' @param start_row int, default 296. 最初进行cgarchfit 所用的数据行数，即in sample data。
  #' @param step_by int. 每隔多久重新fit 模型。
  #' @param multigarchfit_list 保存multigarchfit 的list。
  #' 其fit 的日期必须与copula 模型保持一致，且name 为fit 所在的行数。
  #' @return xts 对象。返回的是每次fit 后向前filter step_by 期的方差-协方差矩阵rcov

  # 数据的总行数，以及需要refit 的行index 数。
  total_rows <- nrow(data)
  fit_time <- seq.int(
    from = start_row,
    to = total_rows, by = step_by
  )

  cls <- parallel::makeForkCluster(parallel::detectCores())
  doParallel::registerDoParallel(cls)
  rolling_rcov <- foreach(t = fit_time, .combine = "rbind") %dopar% {
    # rolling_rcov <- for (t in fit_time) {
    # 1. 解析fit 过的multigarchfit 对象，fit 一个garch-copula 模型

    # 如果t 不在multigarchfit_list 里或者multigarchfit_list 为NULL，都会取出NULL
    fitted_multigarchfit <- multigarchfit_list[[as.character(t)]]

    # 依此使用不同的scale 尝试五次
    for (try_time in 1:5) {
      current_fit <- multisol_cgarchfit(
        spec = pure_cgarch_spec, data = data[1:t, ],
        fit = fitted_multigarchfit,
        fit.control = list(scale = 10 * try_time)
      )

      # 如果拟合返回的不是cGARCHfit 对象，则不使用multigarchfit 对象再来一次
      if (!is(current_fit, "cGARCHfit")) {
        current_fit <- multisol_cgarchfit(
          spec = pure_cgarch_spec, data = data[1:t, ],
          fit = NULL, fit.control = list(scale = 10 * try_time)
        )
      }

      # 如果current_fit 对象类型为cGARCHfit，停止循环即可
      if (is(current_fit, "cGARCHfit")) {
        break
      }
    }

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
    forcasted_cov <- tail(
      rmgarch::rcov(current_filter, output = "matrix"),
      forcast_t
    )
    return(forcasted_cov)
  }
  parallel::stopCluster(cls)

  return(rolling_rcov)
}


rolling_cop_rcov_main <- function() {
  #' @title rolling 计算copula 模型rcov 的主函数
  #' @details 将会分别计算一个t-cop_dcc, norm-cop_dcc, t-cop_static, norm-cop_static
  #' 所对应的rcov，最终保存在一个list 中。

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

  # 1. 读取必要数据
  facs_xts <- read_fac_xts(data_freq = data_freq)
  in_sample_end <- in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS[data_freq])
  multigarchfit_list <- read_rolling_multigarchfit(data_freq = data_freq)

  # 2. 指定cGARCHspec 部分
  arma_order_for_roll <- matrix(rep(3, 10), nrow = 2)
  multi_garch_spec <- all_facs_multigarch(
    arma_order_df = arma_order_for_roll,
    fit = FALSE
  )

  # 3. 正式开始滚动fit
  # 1) t_dcc
  print("Rolling rcov for t-cop dcc.")
  t_dcc_cgarch_spec <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    copula_type = "mvt", is_dcc = TRUE, fit = FALSE
  )
  t_dcc_rcov <- rolling_cgarch_rcov(
    data = facs_xts,
    pure_cgarch_spec = t_dcc_cgarch_spec,
    start_row = in_sample_end,
    step_by = ROLLING_STEP,
    multigarchfit_list = multigarchfit_list
  )

  # 2) norm_dcc
  print("Rolling rcov for norm-cop dcc.")
  norm_dcc_cgarch_spec <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    copula_type = "mvnorm", is_dcc = TRUE, fit = FALSE
  )
  norm_dcc_rcov <- rolling_cgarch_rcov(
    data = facs_xts,
    pure_cgarch_spec = norm_dcc_cgarch_spec,
    start_row = in_sample_end,
    step_by = ROLLING_STEP,
    multigarchfit_list = multigarchfit_list
  )

  # 3) t_static
  print("Rolling rcov for t-cop static.")
  t_static_cgarch_spec <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    copula_type = "mvt", is_dcc = FALSE, fit = FALSE
  )
  t_static_rcov <- rolling_cgarch_rcov(
    data = facs_xts,
    pure_cgarch_spec = t_static_cgarch_spec,
    start_row = in_sample_end,
    step_by = ROLLING_STEP,
    multigarchfit_list = multigarchfit_list
  )

  # 4) norm_static
  print("Rolling rcov for norm-cop static.")
  norm_static_cgarch_spec <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    copula_type = "mvnorm", is_dcc = FALSE, fit = FALSE
  )
  norm_static_rcov <- rolling_cgarch_rcov(
    data = facs_xts,
    pure_cgarch_spec = norm_static_cgarch_spec,
    start_row = in_sample_end,
    step_by = ROLLING_STEP,
    multigarchfit_list = multigarchfit_list
  )

  all_cop_rcov_list <- list(
    t_dcc = t_dcc_rcov, norm_dcc = norm_dcc_rcov,
    t_static = t_static_rcov, norm_static = norm_static_rcov
  )

  saveRDS(all_cop_rcov_list, opts[["out_put"]])
}


if (sys.nframe() == 0) {
  rolling_cop_rcov_main()
}
