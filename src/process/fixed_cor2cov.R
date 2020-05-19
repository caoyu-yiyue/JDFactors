#!/usr/bin/env Rscript

#############################################################################
# 根据multiGARCHfit list fitler 得出的每日sigma，和额外的样本内或外的相关系数
# 矩阵cor matrix，转换计算成“方差-协方差矩阵”cov matrix
# 独立运行脚本时需要传入一个参数，即生成文件的保存路径。
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(rugarch)
  library(MBESS)
})

source("src/data/read_data.R")
source("src/config.R")


in_or_out_sample_cor <- function(data, in_sample_years) {
  #' @title 根据in sample 的年份，对数据的in sample 和out sample 部分分别计算一个相关系数矩阵cor
  #'
  #' @param data xts 对象。需要处理的数据
  #' @param in_sample_years int. 样本内的年数
  #' @return list of 2 cor matrix. 包含两个cor matrix 的list。
  #' names 为c("fixed_cor.IN_SAM", "fixed_cor.OUT_SAM")

  # 先找到in sample 的最后一行
  in_sample_end_row <- in_sample_yearend_row(
    data = data,
    in_sample_year = in_sample_years
  )
  cors_list <- list(
    fixed_cor.IN_SAM = cor(data[1:in_sample_end_row]),
    fixed_cor.OUT_SAM = cor(data[(in_sample_end_row + 1):nrow(data)])
  )
  return(cors_list)
}


rolling_sigma_forcast <- function(multi_garch_fit_list, data, step_by) {
  #' @title 对multiGARCHfit list 滚动filter 预测sigma
  #'
  #' @param multi_garch_fit_list 保存有multiGARCHfit 对象的list，names 为fit 时所用的数据的行数
  #' @param data 进行rolling filter 的整个数据
  #' @param step_by 每隔多少期（即多少行）进行了一次refit
  #' @return xts 对象。是整个样本外时期中filter 出来的sigma。

  total_rows <- nrow(data)
  fit_times <- as.numeric(names(multi_garch_fit_list))

  sigma_result_list <- list()
  for (fit_t in fit_times) {
    # 当前的multigarch fit 对象
    multifit_obj <- multi_garch_fit_list[[as.character(fit_t)]]

    # 当总行数与t 相比大于等于step_by(12)，则filter data 为到t 往后step_by(12)期
    # 预测期为step_by(12)期；否则，filter data 为整个facs_xts，预测期为所有数据行 - t
    if (total_rows - fit_t >= step_by) {
      filter_data <- data[1:(fit_t + step_by), ]
      forcast_t <- step_by
    } else {
      filter_data <- data
      forcast_t <- total_rows - fit_t
    }

    # 进行filter
    current_filt <- multifilter(
      multifitORspec = multifit_obj, data = filter_data,
      n.old = fit_t
    )

    # 输出sigma，然后追加到list 中
    forcasted_sigma <- tail(sigma(current_filt), forcast_t)
    sigma_result_list[[as.character(fit_t)]] <- forcasted_sigma
  }
  result_sigmas <- do.call(rbind, sigma_result_list)

  # 找到预测部分的index 并设置，同时找到原数据的colname 并设置
  sigma_idx <- index(data)[(fit_times[1] + 1):total_rows]
  sigmas_xts <- as.xts(result_sigmas, order.by = sigma_idx)
  colnames(sigmas_xts) <- colnames(data)

  # 返回滚动预测的sigmas
  return(sigmas_xts)
}


cors2covs <- function(cor_mat, sigma_xts) {
  #' @title 根据固定的相关系数矩阵cor_mat，以及每天的sd(这里是sigma_xts)，计算出每天展平的方差-协方差向量
  #' @param cor_mat 一个固定的相关系数矩阵
  #' @param sigma_xts 因子每天的预测sd (即garch 模型预测的sigma)的xts 对象
  #' @details 函数返回的是一个xts 对象，index 与sigma_xts 相同，方差-协方差的顺序为：
  #' 方差-协方差矩阵按行自左向右排列，第二行从该行第二项向右排列，以此类推
  #' 这是使用lower.tri 得到的对称矩阵的向量
  #' @return xts 对象。保存了每日的方差-协方差向量，排列方式如details 所述

  # 定义针对一行sigma 进行cor2cov 变换的函数
  .single_row_cor2cov <- function(sigma_vec, cor_mat) {
    cov_mat <- MBESS::cor2cov(cor.mat = cor_mat, sd = sigma_vec)
    flat_cov_vec <- cov_mat[lower.tri(cov_mat, diag = TRUE)]
    return(flat_cov_vec)
  }

  # 将上面的函数应用到sigma_xts 的每一行上
  covs_mat <- apply(sigma_xts, 1, FUN = .single_row_cor2cov, cor_mat = cor_mat)
  covs_xts <- as.xts(t(covs_mat))

  # 对列名进行组合，加到展开的covs vector 上
  # 算法来自rmgarch 包中未export 的函数make.cov.names
  fac_names <- colnames(sigma_xts)
  n <- length(fac_names)
  names_mat <- matrix(NA, n, n) # 创建一个空matrix，用来放名字组合
  diag(names_mat) <- fac_names
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        next
      }
      names_mat[i, j] <- paste0(fac_names[i], ":", fac_names[j])
    }
  }
  flat_names <- names_mat[lower.tri(names_mat, diag = TRUE)]

  # 命名并输出
  colnames(covs_xts) <- flat_names
  return(covs_xts)
}


fix_cor2cov_main <- function() {
  #' @title 根据multiGARCHfit list 和额外的固定cor matrix 计算每日cov matrix 的主函数
  #' @return NULL. 但将会生成一个list，包含使用样本内和样本外cor matrix 计算所得的
  #' 每日cov matrix 的两个xts 对象。list name 为c("fixed_cor.IN_SAM", "fixed_cor.OUT_SAM")

  # 读取因子数据和multiGARCHfit 对象的list
  facs_xts <- read_fac_xts()
  multigarch_list <- read_rolling_multigarchfit()

  # 计算样本内与样本外的cor 即相关系数矩阵
  cors_list <- in_or_out_sample_cor(facs_xts, IN_SAMPLE_YEARS)
  # 滚动fitler 预测sigma
  forcasted_sigma_xts <- rolling_sigma_forcast(
    multi_garch_fit_list = multigarch_list,
    data = facs_xts,
    step_by = ROLLING_STEP
  )

  # 对样本内与外的cor list 应用函数cors2covs，转换成flat cov xts 对象
  flat_covs_list <- lapply(cors_list,
    FUN = cors2covs,
    sigma_xts = forcasted_sigma_xts
  )

  # 保存到输入的路径当中。
  cmd_args <- commandArgs(trailingOnly = TRUE)
  saveRDS(flat_covs_list, file = cmd_args[[1]])
}


if (sys.nframe() == 0) {
  fix_cor2cov_main()
}
