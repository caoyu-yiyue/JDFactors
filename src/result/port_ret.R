#!/usr/bin/env Rscript

#############################################################################
# 根据计算不同sigma 下不同策略计算所得的权重，得出每日的投资组合收益率序列
#############################################################################

suppressPackageStartupMessages({
  library(xts)
})

source("src/data/read_data.R")


.rets_for_single_strategy <- function(fac_weights, facs_rets) {
  #' @title 对每一个策略计算每日收益率
 
  weighted_facs <- facs_rets * fac_weights
  port_ret <- xts::as.xts(rowSums(weighted_facs),
    order.by = index(weighted_facs)
  )
  return(port_ret)
}


.rets_for_single_sigma <- function(weights_list_in_sigma, facs_ret) {
  #' @title 对每个sigma 下的策略list，计算出一个每日组合收益xts，每列为策略名
 
  strategies_ret_list <- lapply(weights_list_in_sigma, .rets_for_single_strategy, facs_ret)
  strategies_rets_xts <- do.call("cbind", strategies_ret_list)
  colnames(strategies_rets_xts) <- names(weights_list_in_sigma)
  return(strategies_rets_xts)
}


rets_for_all_sigma <- function(opt_weights, facs_ret) {
  #' @title 对opt_weights 二级list(sigma -> strategy 两级)计算每日组合收益序列

  port_rets_list <- lapply(opt_weights, .rets_for_single_sigma, facs_ret = facs_ret)
  return(port_rets_list)
}


invest_result_main <- function() {
  facs_xts <- read_fac_xts()
  opt_weights <- read_opt_weights(which = "all")

  port_rets_list <- rets_for_all_sigma(opt_weights = opt_weights, facs_ret = facs_xts)
  return(port_rets_list)
}
