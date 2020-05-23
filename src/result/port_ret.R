#!/usr/bin/env Rscript

#############################################################################
# 根据计算不同sigma 下不同策略计算所得的权重，得出每日的投资组合收益率序列
# 单独运行脚本时，需要传入一个参数：即生成数据的保存路径。
#############################################################################

suppressPackageStartupMessages({
  library(xts)
})

source("src/data/read_data.R")


.rets_for_single_strategy <- function(fac_weights, facs_rets) {
  #' @title 对每一个策略计算组合的每日收益率
  #'
  #' @param fac_weights xts 对象。各因子的每期权重值
  #' @param facs_rets xts 对象。各因子的每期收益率
  #' @return xts 对象。为通过权重和收益计算所得的组合收益，只有一列，

  weighted_facs <- facs_rets * fac_weights
  port_ret <- xts::as.xts(rowSums(weighted_facs),
    order.by = index(weighted_facs)
  )
  return(port_ret)
}


.rets_for_single_sigma <- function(weights_list_in_sigma, facs_ret) {
  #' @title 对每个sigma 下的策略list，计算出一个每日组合收益xts，每列为策略名
  #'
  #' @param weights_list_in_sigma list. 某个sigma 下的不同策略的权重组成的list
  #' @param facs_ret xts 对象。各因子每期收益率
  #' @return xts 对象。index 与权重相同，col names 为策略名

  strategies_ret_list <- lapply(
    weights_list_in_sigma,
    .rets_for_single_strategy, facs_ret
  )
  strategies_rets_xts <- do.call("cbind", strategies_ret_list)
  colnames(strategies_rets_xts) <- names(weights_list_in_sigma)
  return(strategies_rets_xts)
}


rets_for_all_sigma <- function(opt_weights, facs_ret) {
  #' @title 对opt_weights 二级list(sigma -> strategy 两级)计算每日组合收益序列
  #'
  #' @param opt_weights 保存权重的二级list
  #' @param facs_ret xts 对象。各因子每期收益率
  #' @return list. 每个key 为sigma 值，value 为一个xts 对象，包含了不同策略的组合收益率序列。

  port_rets_list <- lapply(opt_weights,
    .rets_for_single_sigma,
    facs_ret = facs_ret
  )
  return(port_rets_list)
}


invest_result_main <- function() {
  #' @title 计算投资策略的组合收益率结果的主函数
  #' @return NULL，但会将生成的不同sigma 下的各种策略收益率list 保存到传入的第一个参数的位置。

  facs_xts <- read_fac_xts()
  opt_weights <- read_opt_weights(which = "all")

  port_rets_sum1 <- rets_for_all_sigma(
    opt_weights = opt_weights[["sum1"]],
    facs_ret = facs_xts
  )
  port_rets_no_sum1 <- rets_for_all_sigma(
    opt_weights = opt_weights[["no_sum1"]],
    facs_ret = facs_xts
  )

  port_rets_list <- list(sum1 = port_rets_sum1, no_sum1 = port_rets_no_sum1)

  # 保存到输入的路径当中。
  cmd_args <- commandArgs(trailingOnly = TRUE)
  saveRDS(port_rets_list, file = cmd_args[[1]])
}


if (sys.nframe() == 0) {
  invest_result_main()
}
