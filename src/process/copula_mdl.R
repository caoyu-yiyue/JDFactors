#!/usr/bin/env Rscript

#############################################################################
# 拟合copula 模型的函数功能。
# 独立运行脚本时将计算出数据整体的copula 对象, 需要传入参数：保存所以copula 对象的路径
#############################################################################

suppressPackageStartupMessages({
  library(rmgarch)
})

source("src/process/multi_garch_mdl.R")
source("src/data/read_data.R")


# uarch_spec 输入换成multigarch spec 对象
fit_garch_copula <- function(multigarch_spec, copula_type, is_dcc,
                             dcc_order = c(1, 1), asymm = FALSE, fit = TRUE,
                             fac_data = NULL, multigarch_fit = NULL,
                             cluster = NULL) {
  #' @title 基于multigarch_spec 对象和一组数据，计算出一个copula fit 对象
  #'
  #' @param multigarch_spec 一个multigarch spec 对象，将传递给 cgarchspec(uspec)
  #' @param copula_type str, c("mvnorm", "mvt") 的其中之一。copula 模型的类型。
  #' @param is_dcc Bool. 是否使用dcc 即time.varying。
  #' @param asymm Bool, default FALSE. 是否使用非对称的dcc 模型。
  #' @param fit Bool, default TRUE. 是否进行拟合。如果为FALSE 将返回cGARCHspec 对象。
  #' @param fac_data xts 等对象。用于计算模型的数据（xts 等对象）
  #' @param multigarch_fit (optional) uGARCHmultifit 对象，default NULL。
  #' 已经拟合过的multifit 对象。
  #' @param cluster cluster 对象，使用过后需要stopCluster()
  #'
  #' @details 前面的5 个参数均传递给rmgarch::cgarchspec(), 最后两个传递给rmgarch::cgarchfit()
  #' @return fit == TURE 时，cGARCHfit 对象。根据上述的参数拟合完成的cGARCHfit 对象。
  #' fit == FALSE 时，返回cGARCHspec 对象，即不对模型进行拟合而时返回spec。

  # 指定coplua spec
  cop_spec <- cgarchspec(
    uspec = multigarch_spec,
    asymmetric = asymm,
    distribution.model = list(
      copula = copula_type, method = "ML",
      time.varying = is_dcc, transformation = "parametric"
    ),
    dccOrder = dcc_order
  )

  if (!fit) {
    return(cop_spec)
  } else {
    # 进行copula fit
    cop_fit <- cgarchfit(cop_spec,
      data = fac_data, cluster = cluster,
      fit = multigarch_fit
    )

    return(cop_fit)
  }
}


copula_all_main <- function() {
  #' @title 计算整体数据多种copula 对象的主函数
  #' @return 集合了所有copula 的list

  multigarch_fit <- read_multi_garch_fit()
  best_arma_order <- read_best_arma_order(which = "adjusted")
  facs_xts <- read_fac_xts()

  multi_garch_spec <- all_facs_multigarch(
    arma_order_df = best_arma_order,
    fit = FALSE
  )
  cls <- parallel::makeForkCluster(4)

  # 计算每个copula
  static_norm_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvnorm",
    is_dcc = FALSE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  static_t_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvt",
    is_dcc = FALSE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  dcc_norm_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvnorm",
    is_dcc = TRUE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  dcc_t_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvt",
    is_dcc = TRUE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  adcc_norm_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvnorm",
    is_dcc = TRUE,
    asymm = TRUE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  adcc_t_cop <- fit_garch_copula(
    multigarch_spec = multi_garch_spec,
    fac_data = facs_xts,
    copula_type = "mvt",
    is_dcc = TRUE,
    asymm = TRUE,
    multigarch_fit = multigarch_fit,
    cluster = cls
  )
  parallel::stopCluster(cls)

  # 将所有的copula 集合为list
  all_cops <- list(
    static_norm = static_norm_cop, static_t = static_t_cop,
    dcc_norm = dcc_norm_cop, dcc_t = dcc_t_cop,
    adcc_norm = adcc_norm_cop, adcc_t = adcc_t_cop
  )

  args <- commandArgs(trailingOnly = TRUE)
  saveRDS(all_cops, file = args[[1]])
}


if (sys.nframe() == 0) {
  copula_all_main()
}
