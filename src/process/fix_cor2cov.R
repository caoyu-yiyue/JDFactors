#!/usr/bin/env Rscript

#############################################################################
# 使用样本内或样本外的相关系数矩阵，结合garch 模型计算出的sigma，得到固定相关
# 系数矩阵的方差-协方差矩阵
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(rugarch)
  library(foreach)
})

source("src/data/read_data.R")
source("src/config.R")


rolling_sigma_forcast <- function(multi_garch_fit_list, data, step_by) {
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


fix_cor2cov_main <- function() {
  facs_xts <- read_fac_xts()
  multigarch_list <- read_rolling_multigarchfit()

  forcasted_sigma_xts <- rolling_sigma_forcast(
    multi_garch_fit_list = multigarch_list,
    data = facs_xts,
    step_by = ROLLING_STEP
  )
}

fix_cor2cov_main()
