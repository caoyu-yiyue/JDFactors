#!/usr/bin/env Rscript

#############################################################################
# 对每个因子进行多种arma order 的NGARCH 拟合，选出最小的BIC 的arma order 并保存
# 单独执行该文件时，需要传入的参数为：
#   --data_freq 数据频率，Week, Month 或Day，即哪个频率的因子数据
#   -o output file path
#############################################################################


suppressPackageStartupMessages({
  library(rugarch)
  library(future.apply)
  library(optparse)
})
source("src/data/read_data.R")


single_arma_garch_ais <- function(arma_order, data_vector,
                                  distribution = "sstd") {
  #' @title 针对一个arma order，对一个单维数据计算出BIC
  #'
  #' @param arma_order 两个整数组成的numeric，指定arma order
  #' @param data_vector 一个单维数据vector，计算garch 模型的数据
  #' @param distribution defalut "sstd". garch 模型的假定分布
  #' @return 返回一个整数，即garch 模型计算所得的BIC 值
  #' @details 如果arma_order 的和为0，函数将返回一个Inf，以示该情况不合法

  # 如果arma order 和为0，则返回Inf。即舍弃此种情况。
  if (sum(arma_order) == 0) {
    return(Inf)
  }

  # 指定garch 模型
  var_mdl <- list(model = "fGARCH", garchOrder = c(1, 1), submodel = "NGARCH")
  mean_mdl <- list(armaOrder = arma_order)
  garch_spec <- ugarchspec(
    variance.model = var_mdl, mean.model = mean_mdl,
    distribution.model = distribution
  )
  # 拟合并取出BIC 值
  fit <- ugarchfit(spec = garch_spec, data = data_vector, solver = "hybrid")
  bic <- infocriteria(fit)["Bayes", ]
  return(bic)
}


single_fac_best_arma_order <- function(data_vector, arma_order_pairs_df,
                                       arch_dist = "sstd") {
  #' @title 对一个单维数据，输入一组arma order(每行一组)，返回BIC 最大的arma order 值
  #'
  #' @param data_vector 用于计算garch 模型的单维数据
  #' @param arma_order_pairs_df 每行一组arma order 的data.frame
  #' @param arch_dist arch 模型的假定分布
  #' @return 一个两个整数的numeric vector，为最佳的arma order
  #' @details 该函数使用future.apply::future_apply 进行并行计算，需要在函数外通过
  #' future.apply::plan() 函数指定并行运算方式和workers 数量。

  # 对每个arma order 循环，计算输入单列数据的最佳BIC 并返回
  all_bic <- future_apply(arma_order_pairs_df, 1, single_arma_garch_ais,
    data_vector = data_vector, distribution = arch_dist
  )
  best_arma_order <- as.numeric(arma_order_pairs_df[which.min(all_bic), ])
  return(best_arma_order)
}


best_arma_mian <- function(fac_data_freq, save_path) {
  #' @title 计算所有因子最佳arma order 的主函数
  #' @param fac_data_freq 哪个频率的因子数据 c("Week", "Day", "Month")
  #' @param save_path 保存最佳arma order 数据的文件地址（完整文件名）
  #' @return NULL
  plan(multiprocess, workers = 4)
  fac_xts <- read_fac_xts(data_freq = fac_data_freq)
  arma_orders <- expand.grid(ar = 0:10, ma = 0:10)
  best_arma_orders <- apply(fac_xts, 2, single_fac_best_arma_order,
    arma_order_pairs_df = arma_orders
  )
  saveRDS(best_arma_orders, file = save_path)
}


if (!interactive()) {
  # 读取命令行输入
  option_list <- list(
    make_option(
      opt_str = "--data_freq", type = "character",
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

  # 必须提供文件存储路径
  if (is.null(opts[["out_put"]])) {
    print_help(opt_parser)
    stop("Must input the out put file path")
  }

  # 运行主函数
  best_arma_mian(
    fac_data_freq = opts[["data_freq"]],
    save_path = opts[["out_put"]]
  )
}
