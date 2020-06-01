#!/usr/bin/env Rscript

#############################################################################
# 根据multiGARCHfit list fitler 得出的每日sigma，和额外的样本内或外的相关系数
# 矩阵cor matrix，转换计算成“方差-协方差矩阵”cov matrix。
# 此外，再计算一个使用样本内的cov 在样本外时期重复的benchmark cov 矩阵。
# 独立运行脚本时需要传入一个参数，即生成文件的保存路径。
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(rugarch)
  library(MBESS)
  library(optparse)
})

source("src/data/read_data.R")
source("src/config.R")


in_or_out_sample_cor <- function(data, in_sample_end_row) {
  #' @title 根据in sample 的年份，对数据的in sample 和out sample 部分分别计算一个相关系数矩阵cor
  #'
  #' @param data xts 对象。需要处理的数据
  #' @param in_sample_years int. 样本内的年数
  #' @return list of 2 cor matrix. 包含两个cor matrix 的list。
  #' names 为c("fixed_cor.IN_SAM", "fixed_cor.OUT_SAM")

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


.make_flat_cov_names <- function(fac_names) {
  #' @title 对输入的列名进行组合，然后按照cov 矩阵lower.tri 的方向展开成向量。
  #' @param fac_names 因子名，即需要进行组合的列名
  #' @details 算法来自rmgarch 包中未export 的函数make.cov.names
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
  return(flat_names)
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

  # 将上面的函数应用到sigma_xts 的每一行上，结果需要转置，并且重新设置xts 对象
  # 直接对apply 的结果转置再生成xts 将导致index 带上时区，以后很难对齐数据
  covs_mat <- apply(sigma_xts, 1, FUN = .single_row_cor2cov, cor_mat = cor_mat)
  covs_xts <- as.xts(t(covs_mat), order.by = as.Date(colnames(covs_mat)))

  fac_names <- colnames(sigma_xts)
  flat_names <- .make_flat_cov_names(fac_names = fac_names)

  # 命名并输出
  colnames(covs_xts) <- flat_names
  return(covs_xts)
}


static_in_sam_cov <- function(data, in_sample_end_row) {
  #' @title 计算样本期内的cov，展平并填充到样本外的每一期，作为benchmark 使用
  #'
  #' @param data 计算cov 的数据
  #' @param in_sample_end_row in sample 部分的行数
  #' @return xts 对象。index 为样本外时间，每一行都是样本内的展平cov。

  total_row <- nrow(data)
  in_sample_data <- data[1:in_sample_end_row, ]
  out_sam_index <- index(data)[(in_sample_end_row + 1):total_row]

  # 计算 in_sample 时期的cov，并展平
  in_sam_cov <- cov(in_sample_data)
  flat_cov <- in_sam_cov[lower.tri(in_sam_cov, diag = TRUE)]

  # 将一行重复out_of_sample 行数一样的次数，组成一个matrix
  flat_cov_mat <- do.call("rbind", replicate(
    n = length(out_sam_index),
    flat_cov, simplify = FALSE
  ))

  # 给matrix 的列重命名，转换为xts 同时加入index
  colnames(flat_cov_mat) <- .make_flat_cov_names(colnames(data))
  static_cov_xts <- as.xts(flat_cov_mat, order.by = out_sam_index)

  return(static_cov_xts)
}


fix_cor2cov_main <- function() {
  #' @title 根据multiGARCHfit list 和额外的固定cor matrix 计算每日cov matrix 的主函数
  #' @return NULL. 但将会生成一个list，包含使用样本内和样本外cor matrix 计算所得的
  #' 每日cov matrix 的两个xts 对象。
  #' list name 为c("fixed_cor.IN_SAM", "fixed_cor.OUT_SAM")

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

  # 读取因子数据和multiGARCHfit 对象的list
  facs_xts <- read_fac_xts(data_freq = data_freq)
  multigarch_list <- read_rolling_multigarchfit(data_freq = data_freq)
  in_sample_end_row <- in_sample_yearend_row(
    data = facs_xts,
    in_sample_year = IN_SAMPLE_YEARS[data_freq]
  )

  # 计算样本内与样本外的cor 即相关系数矩阵
  cors_list <- in_or_out_sample_cor(facs_xts, in_sample_end_row)
  # 滚动fitler 预测sigma
  forcasted_sigma_xts <- rolling_sigma_forcast(
    multi_garch_fit_list = multigarch_list,
    data = facs_xts,
    step_by = ROLLING_STEP[data_freq]
  )

  # 对样本内与外的cor list 应用函数cors2covs，转换成flat cov xts 对象
  flat_covs_list <- lapply(cors_list,
    FUN = cors2covs,
    sigma_xts = forcasted_sigma_xts
  )

  # 计算静态的in sample cov，并加入到结果list 中
  static_cov_xts <- static_in_sam_cov(
    data = facs_xts,
    in_sample_end_row = in_sample_end_row
  )
  flat_covs_list[["static_benchmark"]] <- static_cov_xts

  # 保存到输入的路径当中。
  saveRDS(flat_covs_list, file = opts[["out_put"]])
}


if (sys.nframe() == 0) {
  fix_cor2cov_main()
}
