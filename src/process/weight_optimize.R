#!/usr/bin/env Rscript

#############################################################################
# 使用rolling mean 和rolling rcov，解最优化求得最优权重的模块儿
# 传入一个参数：保存结果list 的路径
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(quadprog)
})

source("src/data/read_data.R")


.vec2var_mat <- function(var_vec, n_fac) {
  #' @title 从一个序列的var-cvar 构造出一个var-cvar 矩阵
  #'
  #' @param var_vec vector. var_cvar 序列。顺序为：按协方差矩阵的行排列，第一行全部包括
  #' 第二行包括2 到最后一项，以此类推。
  #' @param n_fac 因子的个数，即协方差矩阵的行（列）数
  #' @details var_vec 为按行组成var-cvar 矩阵的右上角。
  #' @return matrix 方差-协方差矩阵

  m1 <- matrix(NA, n_fac, n_fac)
  m1[lower.tri(m1, diag = TRUE)] <- var_vec
  m2 <- t(m1)
  m2[lower.tri(m2, diag = TRUE)] <- var_vec
  return(m2)
}


.opt_single_day <- function(mean_var_row, risk_pref_coef,
                            amat, bvec, meq, n_fac = 5) {
  #' @title 用于对每一行数据进行最优化的函数
  #'
  #' @param mean_var_row numeric 一行数据，包括了前n 行的mean 和
  #' 后面n * (n + 1) / 2 的var-cvar 数列
  #' @param risk_pref_coef scalar 风险偏好系数gamma
  #' @param amat matrix 传入quadprog::solve.QP 的Amat
  #' @param bvec vecotr 传入quaprog::solve.QP 的bvec
  #' @param meq scalar 传入quaprog::solve.QP 的meq
  #' @param n_fac scalar, default 5. 因子的数量。
  #' @return numeric vector. 最优权重。

  # 分别取出mean 和var-cvar 部分，并将cvar 部分转化为matrix
  mean_vec <- mean_var_row[1:n_fac]
  var_cvar_vec <- mean_var_row[(n_fac + 1):length(mean_var_row)]
  var_cvar_mat <- .vec2var_mat(var_vec = var_cvar_vec, n_fac = n_fac)
  # 将var-cvar matrix 乘以系数
  opt_mat <- risk_pref_coef * var_cvar_mat

  # 解最优化
  sol_result <- solve.QP(
    Dmat = opt_mat, dvec = mean_vec,
    Amat = amat, bvec = bvec, meq = meq
  )
  return(sol_result[["solution"]])
}


rolling_opt <- function(data, n_fac = 5, sum_1 = TRUE) {
  #' @title 针对一个mean 和var-cvar 合并在一起的数据框，滚动求出最优权重
  #' @param data 前n_fac 列为mean，后面的n * (n + 1) / 2 列为var-cvar 数列的数据
  #' 其中，mean 列的列名将成为权重列的列名。
  #' @param n_fac 因子（mean 数据）的个数
  #' @param sum_1 是否带有权重和为1 的限制条件。
  #' @details 将对传入的数据进行二次最优化问题求解。必然带上每个权重大于等于0 的条件，
  #' 同时通过sum_1 确定是否加上权重和为1 的条件。
  #' @return 最优化后所得的最优权重值。

  # 对于因子和是否为1 进行判断，给出不同的限制条件
  if (sum_1) {
    amat <- cbind(
      matrix(rep.int(1, n_fac), nrow = n_fac), # 权重和为1
      diag(n_fac) # 每个权重大于等于0
    )
    bvec <- c(1, rep.int(0, n_fac))
    meq <- 1
  } else {
    amat <- diag(n_fac)
    bvec <- rep.int(0, n_fac)
    meq <- 0
  }

  # 使用上面的限制条件，对每一行解最优化
  weights_matrix <- apply(data,
    MARGIN = 1,
    FUN = .opt_single_day, risk_pref_coef = 7,
    amat = amat, bvec = bvec, meq = meq, n_fac = 5
  )

  # 整理apply 传回的数据的格式，转置、变xts、加名字
  trans_xts_obj <- as.xts(t(weights_matrix))
  names(trans_xts_obj) <- names(data)[1:n_fac]

  return(trans_xts_obj)
}


roll_opt_main <- function() {
  #' @title 进行权重最优化的主函数

  # 读取mean 数据
  rolling_mean <- read_rolling_mean()
  n_fac <- ncol(rolling_mean)

  # 读取保存所有rcov 的list
  all_rcovs <- read_rolling_cop_rcov(which = "all")
  rcov_names <- names(all_rcovs)

  # 对每个rcov 循环，与mean 合并并最优化
  opt_weights <- list()
  for (name in rcov_names) {
    rcov_xts <- read_rolling_cop_rcov(which = name)
    mean_rcov_merged <- cbind(rolling_mean, rcov_xts)
    opt_weight <- rolling_opt(data = mean_rcov_merged, n_fac = n_fac, sum_1 = TRUE)
    opt_weights[[name]] <- opt_weight
  }

  # 保存到输入的路径当中。
  cmd_args <- commandArgs(trailingOnly = TRUE)
  saveRDS(opt_weights, file = cmd_args[[1]])
}


if (sys.nframe() == 0) {
  test_list <- roll_opt_main()
}
