#!/usr/bin/env Rscript

#############################################################################
# 滚动计算每一期哪个arma order 有最小的bic
#############################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(rugarch)
  library(future.apply)
})

source("src/config.R")
source("src/data/read_data.R")


ugarchfit_retry <- function(spec, data, use_default_par = TRUE,
                            retry_time = 3) {
  #' @title 一个重试ugarchfit 函数的wrapper 函数。重试(retry_time * 2 + 1) 次
  #' 返回第一个合法结果，或返回最后一个结果
  #' @param spec uGARCHspec，指定GARCH 模型的设定
  #' @param data 用于拟合garch 的数据。这里必须是单维数据（用于ugarchfit 函数)。
  #' @param use_default_par Bool, 默认TRUE. 是否使用默认参数先fit 一次。
  #' @param retry_time 重试次数，默认3。这里首先使用默认参数的"hybrid" solver 计算一次
  #' 然后，使用不同的scale 参数重试retry_time 次。最后再使用"lbfgs" solver 重试
  #' retry_time 次，共(retry_time * 2 + 1)次
  #' @return uGARCHfit 对象。第一个成功通过检验的uGARCHfit 对象；或重试最后一次的uGARCHfit

  # 如果指定使用默认参数拟合，首先使用默认参数进行一次拟合
  if (use_default_par) {
    ugfit <- ugarchfit(spec = spec, data = data, solver = "hybrid")
    if (ugfit@fit$convergence == 0 & "cvar" %in% names(ugfit@fit)) {
      # 通过验证，返回结果
      return(ugfit)
    }
  }

  # 开始重试
  for (try_time in 1:retry_time) {
    ugfit <- ugarchfit(
      spec = spec, data = data, solver = "hybrid",
      fit.control = list(scale = 10**try_time)
    )

    if (ugfit@fit$convergence == 0 & "cvar" %in% names(ugfit@fit)) {
      # 通过验证，返回结果
      return(ugfit)
    }

    # 如果到了最后一次还没有成功，会运行到这里，尝试solver：lbfgs 五次
    for (lbfgs_try_time in 1:retry_time) {
      ugfit <- ugarchfit(
        spec = spec, data = data, solver = "lbfgs",
        fit.control = list(scale = 10**try_time)
      )
      if (ugfit@fit$convergence == 0 & "cvar" %in% names(ugfit@fit)) {
        # 通过验证，返回结果
        return(ugfit)
      }
    }
  }

  # 如果以上的10 次尝试都失败，会走到这里。返回最后一次的ugfit
  return(ugfit)
}


single_arma_garch_bis <- function(info_row, full_data,
                                  distribution = "sstd") {
  #' @title 应用在一个arma order combine 数据框每一行上的函数。
  #' 对每一行拟合一个uGARCHfit，返回BIC 值和conv 检验的flag
  #' @param info_row numeric. 一个保存了数据位置和arma order 的row。共四个元素。
  #' 分别为数据行坐标、数据列坐标、ar order、ma order
  #' @param full_data data.frame/xts. 一个完整的数据框。
  #' 使用info_row 的前两个元素找到本次fit 的单维数据
  #' @param distribution str, default "sstd". garch 模型的假定分布，传给ugarchspec。
  #' @return two element numeric. 第一个元素为ugarchfit 的bic 值；
  #' 第二个为conv test 的flag 值（由bool 转为int）

  # 解析arma order 和用于本次拟合的数据data_vector
  arma_order <- info_row[3:4]
  data_vector <- full_data[1:info_row["t"], info_row["fac_idx"]]

  # 指定garch 模型
  var_mdl <- list(model = "fGARCH", garchOrder = c(1, 1), submodel = "NGARCH")
  mean_mdl <- list(armaOrder = arma_order)
  garch_spec <- ugarchspec(
    variance.model = var_mdl, mean.model = mean_mdl,
    distribution.model = distribution
  )
  # 拟合garch 模型
  fit <- ugarchfit_retry(spec = garch_spec, data = data_vector, retry_time = 3)

  # 进行conv test
  conv <- !as.logical(fit@fit$convergence)
  have_cvar_flag <- "cvar" %in% names(fit@fit)
  conv_flag <- as.numeric(conv & have_cvar_flag)

  bic <- tryCatch(infocriteria(fit)["Bayes", ], error = function(cond) NA)
  return(c(bic, conv_flag))
}


extract_best_armas <- function(bic_and_conv, arma_order_combine) {
  #' @title 对每天的数据，找到通过conv 检验后的最小BIC 的arma order
  #' @param bic_and_conv 一个两列的数据框。分别为bic 值和conv flag
  #' @param arma_order_combine 起初计算bic 和conv 的arma 数据框。四列分别为
  #' 时间（行号）、列号、ar order、ma order
  #' @return dplyr::tbl 对象。保存每个拟合时点、每个因子的最佳arma order。
  #' 一个(2 + n_fac) 列的数据框，第一列为原数据时间（行号）
  #' 第二列为ar/ma 名字、后n_fac 列为对应的因子序号。

  # 先对bic 和conv 命名列，然后与原arma 数据合并
  colnames(bic_and_conv) <- c("bic", "conv_flag")
  merged_result <- cbind(arma_order_combine, bic_and_conv)

  # 留下非空的bic 以及conv_flag 解决的行
  result_tbl <- as.tbl(merged_result)
  filt_tbl <- result_tbl %>%
    filter(!is.na(bic) & conv_flag == 1) %>%
    group_by(t, fac_idx) %>%
    slice(which.min(bic)) %>%
    ungroup()

  # 转换数据格式
  long_tbl <- tidyr::pivot_longer(
    data = filt_tbl[, 1:4],
    cols = c("ar", "ma")
  ) %>%
    tidyr::pivot_wider(names_from = "fac_idx")
  return(long_tbl)
}


current_arma_order <- function(best_armas_df, t) {
  #' @title 根据计算所得的所有时间best arma df，找到t 时刻的最佳arma order
  #' @param best_armas_df tbl(df) 保存所有时刻最佳arma order 的tbl(df) 对象
  #' @param t numeric. 要找的时间（行号）
  #' @return matrix. t 时刻的最佳arma order。如果该时刻没有计算过，会返回一个
  #' 0 行的matrix。

  current_armas <- best_armas_df[best_armas_df$t == t, 3:ncol(best_armas_df)]
  current_armas_mat <- as.matrix(current_armas)
  return(current_armas_mat)
}


rolling_best_arma_main <- function() {
  option_list <- list(
    make_option(
      opt_str = c("-f", "--data_freq"), type = "character",
      default = "Week", help = "Which data freq of factors?",
      metavar = "character"
    ),
    make_option(c("-o", "--out_put"),
      type = "character",
      default = NULL, help = "Path to save best arma order data.",
      metavar = "character"
    )
  )
  opt_parser <- OptionParser(option_list = option_list)
  opts <- parse_args(opt_parser)
  data_freq <- opts[["data_freq"]]

  # 读取数据
  facs_xts <- read_fac_xts(data_freq = data_freq)
  # 找到配置参数
  in_sample_end <- if (data_freq == "Month") {
    IN_SAMPLE_YEARS[data_freq]
  } else {
    in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS[data_freq])
  }
  fit_times <- seq.int(in_sample_end, nrow(facs_xts), ROLLING_STEP[data_freq])

  # arma order 的组合数据框
  arma_order_combine <- as.tbl(
    expand.grid(t = fit_times, fac_idx = 1:5, ar = 0:3, ma = 0:3)
  )
  # #排序并去掉ar ma 都是0 的行。
  arma_order_combine <- arma_order_combine %>%
    arrange(t, fac_idx) %>%
    dplyr::filter(ar != 0 | ma != 0) %>%
    as.data.frame()

  # 开始主要计算步骤
  plan(multiprocess, workers = 4)
  bics_wide <- future_apply(arma_order_combine, 1,
    single_arma_garch_bis,
    full_data = facs_xts
  )
  bics_result <- t(bics_wide)
  # 临时保存结果，防止后面出错前面白算。
  best_armas <- extract_best_armas(
    bic_and_conv = bics_result,
    arma_order_combine = arma_order_combine
  )

  rolling_order_list <- list(best_armas = best_armas, all_bics = bics_result)
  saveRDS(rolling_order_list, file = opts[["out_put"]])
}


if (sys.nframe() == 0) {
  rolling_best_arma_main()
}
