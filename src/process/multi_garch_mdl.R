#!/usr/bin/env Rscript

#############################################################################
# 读取中间数据的接口
# 运行脚本时需要输入一个参数：multi garch fit 对象保存的位置
#############################################################################

suppressPackageStartupMessages({
  library(rugarch)
  library(rmgarch)
  library(parallel)
})

source("src/data/read_data.R")


.single_garch_spec <- function(arma_order, garch_order) {
  #' @title 输入一个arma_order(两个整数), 返回一个单变量arch_spec
  #'
  #' @param arma_order arma 阶数，两个整数的numeric vector, 比如c(1, 1)
  #' @return rugarch::uGARCHspec object

  var_mdl <- list(
    model = "fGARCH", garchOrder = garch_order,
    submodel = "NGARCH"
  )
  mean_mdl <- list(armaOrder = arma_order)
  arch_spec <- ugarchspec(
    variance.model = var_mdl, mean.model = mean_mdl,
    distribution.model = "sstd"
  )
  return(arch_spec)
}


all_facs_multigarch <- function(arma_order_df,
                                garch_order_df = matrix(rep(1, 10), nrow = 2),
                                facs_data = NULL, fit = TRUE) {
  #' @title 输入包含每个变量的arma order matrix(df)，对该数据进行多变量garchfit
  #'
  #' @param arma_order_df 保存每个变量的arma order 的matrix(df), 每列为一个fac 的arma order
  #' @param facs_data data.frame, default NULL 需要进行garch fit 的因子数据
  #' @param fit Bool, default TRUE, 指定是否需要进行fit。
  #' @return fit == TURE: rugarch::uGARCHmultifit 即数据每列进行fit 之后的multifit
  #' fit == FALSE: rugarch::uGARCHmultispec 对象，即garch multi spec

  # 对传入的arma_order_df 每列应用.single_garch_spec, 生成一个uGARCHspec list
  garch_speclist <- list()
  for (i in seq_len(ncol(arma_order_df))) {
    arch_spec <- .single_garch_spec(
      arma_order = arma_order_df[, i],
      garch_order = garch_order_df[, i]
    )
    garch_speclist <- append(garch_speclist, arch_spec)
  }
  # garch_speclist <- apply(arma_order_df, 2, .single_garch_spec)
  multi_garch_spec <- multispec(garch_speclist)

  # 如果fit 为FLASE，那么直接返回multi_spec obj
  if (!fit) {
    return(multi_garch_spec)
  }

  # multi fit
  cls4 <- makeCluster(4, type = "FORK")
  multi_fit <- multifit(
    multispec = multi_garch_spec, data = facs_data,
    solver = "nlminb", cluster = cls4
  )
  stopCluster(cls4)
  return(multi_fit)
}


.conver_for_multigarchfit <- function(mulgfit_obj) {
  #' @title 判定一个multigarch_fit 对象中，是否每一个都解开了最优化
  #' @param mulgfit_obj rugarch::multigarchfit object。
  #' @return vector of Bool 标示每个变量是否都conver的flag

  ugfit_objs <- mulgfit_obj@fit
  conver_flags <- sapply(ugfit_objs, function(ugfit) ugfit@fit$convergence)
  bool_conver_flags <- !as.logical(conver_flags)
  return(bool_conver_flags)
}


.cvar_for_multigarchfit <- function(mulgfit_obj) {
  #' @title 判定一个multigarch_fit 对象中，是否每一个ugarchfit@fit 都有slot cvar
  #' @param mulgfit_obj rugarch::multigarchfit object。
  #' @return vector of Bool 标示每个rugarchfit 中是否都有"cvar"的flag

  ugfit_objs <- mulgfit_obj@fit
  have_cvar_flags <- sapply(ugfit_objs, function(ugfit) {
    "cvar" %in% names(ugfit@fit)
  })
  return(as.logical(have_cvar_flags))
}


getspec_multifit <- function(mgfit) {
  #' @title 对一个uGARCHmultifit 对象，找到它的uGARCHmultispec 对象
  #' @param mgfit multiGARCHfit 对象
  #' @return rugarch::uGARCHmultispec 对象

  ugspec_list <- lapply(mgfit@fit, rugarch::getspec)
  result_mspec <- rugarch::multispec(ugspec_list)
  return(result_mspec)
}


multi_garch_fit_main <- function() {
  #' @title 拟合mutliGARCHfit 的主函数
  #'
  #' @return rugarch::uGARCHmultifit 即数据没列进行fit 之后的multifit

  best_arma <- read_best_arma_order(which = "adjusted")
  facs_xts <- read_fac_xts()
  multi_fit_obj <- all_facs_multigarch(
    arma_order_df = best_arma,
    facs_data = facs_xts,
    fit = TRUE
  )

  # 读取命令行输入的唯一参数，保存multi fit obj 文件
  args <- commandArgs(trailingOnly = TRUE)
  saveRDS(multi_fit_obj, file = args[[1]])
}


if (sys.nframe() == 0) {
  multi_garch_fit_main()
}
