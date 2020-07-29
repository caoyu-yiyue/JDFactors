
basic_desc <- function(facs_xts, rf = 0, row_names = NULL) {
  #' @title 对因子数据的基础描述性统计，年化均值、年化方差、夏普比率、
  #' 偏度、峰度。
  #'
  #' @param facs_xts xts 对象，保存因子的历史收益
  #' @param rf xts 对象，保存无风险收益率
  #' @param row_names character. 行名，指定统计量的名称
  #' @return data.frame 即描述性统计df

  library(PerformanceAnalytics)

  ret_df <- table.AnnualizedReturns(facs_xts, Rf = rf)
  min_df <- t(as.data.frame(apply(facs_xts, 2, min)))
  max_df <- t(as.data.frame(apply(facs_xts, 2, max)))
  skew_df <- skewness(facs_xts)
  kurt_df <- kurtosis(facs_xts)
  desc_df <- rbind(ret_df, min_df, max_df, skew_df, kurt_df)
  desc_df <- round(desc_df, 4)
  if (!is.null(row_names)) {
    rownames(desc_df) <- row_names
  }

  return(desc_df)
}


combine_coef_pvalue <- function(ugfit_obj) {
  matcoef <- ugfit_obj@fit$matcoef
  coeffs <- matcoef[, 1]
  pvalues <- matcoef[, 4]

  sig_stars <- gtools::stars.pval(p.value = pvalues)
  pvalues_str <- sprintf("%.4f", pvalues)
  p_with_star <- paste0("(", pvalues_str, sig_stars, ")")

  coef_str <- sprintf("%.4f", coeffs)
  coef_with_pvalue <- paste(coef_str, p_with_star, sep = "\\\n")

  coef_p_df <- as.data.frame(t(as.data.frame(coef_with_pvalue)))
  colnames(coef_p_df) <- names(coeffs)
  rownames(coef_p_df) <- NULL

  # names(coef_with_pvalue) <- names(coeffs)
  return(coef_p_df)
}


weights_LB_for_ugfit <- function(ugfit_obj) {
  #' @title 计算一个uGARCHfit 对象的标准化残差和其平方的的 weighted LB 检验的 p 值
  #'
  #' @param ugfit_obj rugarch::uGARCHfit 对象
  #' @return vector LB 和 LB^2 的 p 值

  std_resd <- rugarch::residuals(ugfit_obj, standardize = TRUE)
  modelinc <- ugfit_obj@model$modelinc

  # 计算不平方的LB 检验
  df_1 <- sum(modelinc[c("ar", "ma")])
  LB_1 <- WeightedPortTest::Weighted.Box.test(std_resd,
    lag = 20,
    type = "Ljung-Box", fitdf = df_1
  )
  p_value_1 <- LB_1[["p.value"]]

  # 计算平方后的LB 检验
  df_2 <- sum(modelinc[c("alpha", "beta")])
  LB_2 <- WeightedPortTest::Weighted.Box.test(std_resd,
    lag = 20,
    type = "Ljung-Box", fitdf = df_2,
    sqrd.res = TRUE
  )
  p_value_2 <- LB_2[["p.value"]]

  p_values <- c(p_value_1, p_value_2)
  names(p_values) <- c("LB", "LB_2")
  return(p_values)
}


.coef_append_pvalue <- function(coefs, pvalues) {
  stars_str <- gtools::stars.pval(pvalues)
  pvalues_str <- sprintf("%.4f", pvalues)
  p_with_star <- paste0("(", pvalues_str, stars_str, ")")

  coef_str <- sprintf("%.4f", coefs)
  coef_with_p <- paste(coef_str, p_with_star, sep = "\\\n")
  return(coef_with_p)
}


cop_param_table <- function(cgfit_obj, fac_names) {
  #' @title 输入一个静态cGARCHfit obj，返回描述所用的表格。
  #'
  #' @param cgfit_obj 一个cGARCHfit 对象
  #' @param fac_names charactor vector 因子的名称
  #' @return 返回一个描述统计表，上方是相关系数参数，下方是自由度、诊断值等

  fac_num <- length(fac_names)
  corr_param_num <- fac_num * (fac_num - 1) / 2
  n_obs <- cgfit_obj@model$modeldata$T # 观测的数量，即拟合时数据的行数
  cop_type <- cgfit_obj@model$modeldesc$distribution # 是norm 还是t copula
  dcc <- cgfit_obj@model$modeldesc$timecopula # 是动态（dcc)还是静态

  # 计算相关系数参数矩阵 corr_param
  cop_param <- coef(cgfit_obj, type = "dcc")
  if (!dcc) {
    corr_param <- copula::p2P(cop_param[1:corr_param_num])
  } else {
    z <- cgfit_obj@mfit$Z
    mean_Q_bar <- (t(z) %*% z) / n_obs
    diag_Q <- diag(diag(mean_Q_bar))
    inverse_diag_Q_sqrt <- sqrt(matrixcalc::matrix.inverse(diag_Q))
    corr_param <- inverse_diag_Q_sqrt %*% mean_Q_bar %*% inverse_diag_Q_sqrt
  }
  corr_param[lower.tri(corr_param, diag = TRUE)] <- NA
  corr_param <- round(corr_param[-fac_num, ], 4) # 最后一行都是NA 不需要

  # 计算其他cop 指标
  if (cop_type == "mvt") {
    shape_param <- sprintf("%.4f", cop_param["mshape"])
    shape_row <- c("$v$", shape_param, rep("", fac_num + 1 - 2))
  }
  if (dcc) {
    # 对于dcc 模型，需要找到a, b, g 和对应的p 值
    cgfit_matcoef <- cgfit_obj@mfit$matcoef
    dcc_matcoef <- cgfit_matcoef[grepl(".+dcc.+", rownames(cgfit_matcoef)), ]
    dcc_param <- .coef_append_pvalue(
      coefs = dcc_matcoef[, 1],
      pvalues = dcc_matcoef[, 4]
    )

    dcc_row <- if (length(dcc_param) == 2) {
      c("a", dcc_param[1], "b", dcc_param[2], rep("", fac_num + 1 - 4))
    } else if (length(dcc_param) == 3) {
      c("a", dcc_param[1], "b", dcc_param[2], "g", dcc_param[3])
    }
  }

  other_cop_params <- if (cop_type == "mvnorm" & !dcc) {
    NULL
  } else if (cop_type == "mvt" & !dcc) {
    shape_row
  } else if (cop_type == "mvnorm" & dcc) {
    dcc_row
  } else if (cop_type == "mvt" & dcc) {
    rbind(shape_row, dcc_row)
  }

  # 诊断指标计算
  llh <- sprintf("%.4f", likelihood(cgfit_obj), 4)
  inf_tests <- rugarch:::.information.test(
    LLH = cgfit_obj@mfit$llh, nObs = n_obs,
    nPars = length(cgfit_obj@mfit$matcoef[, 1])
  )
  BIC <- sprintf("%.4f", inf_tests[["BIC"]], 4)
  diagnosis <- c("LLH", llh, "BIC", BIC, rep("", fac_num + 1 - 4))

  # 排版结果
  result_tbl <- cbind(fac_names[1:(fac_num - 1)], corr_param)
  colnames(result_tbl) <- c("", fac_names)

  result_tbl <- if (is.null(other_cop_params)) {
    rbind(result_tbl, diagnosis)
  } else {
    rbind(result_tbl, other_cop_params, diagnosis)
  }
  rownames(result_tbl) <- NULL
  return(result_tbl)
}


.format_result_table <- function(result_table, annualize_scale) {
  #' @title 对每个gamma 下的结果表格进行需要的格式整理，
  #' 选择需要的策略、对百分数乘100，计算年化的delta U。
  #'
  #' @param result_table matrix/data.frame 单个gamma 下的结果表格。
  #' @param annualize_scale 年化时使用的scale 值。
  #' @return matrix 处理过后的结果表格。

  reranged_tbl <- result_table[c(
    "t_dcc", "norm_dcc", "t_static", "norm_static",
    "sample", "cop_cor.OUT_SAM", "cop_cor.IN_SAM"
  ), ]
  diffed_u <- (reranged_tbl[, "U_n"] - reranged_tbl["cop_cor.IN_SAM", "U_n"])
  delta_n <- (1 + diffed_u)**annualize_scale - 1
  all_cols_df <- cbind(reranged_tbl, delta_n)

  # 每行乘放大系数
  scale_vec <- c(100, 100, 1, 1, 1, 10000, 1, 1, 10000, 100)
  scaled_df <- t(apply(all_cols_df, MARGIN = 1, function(row) scale_vec * row))

  # 保留两位小数并保存为str
  result_df <- format(round(scaled_df, digits = 2), nsmall = 2)
  result_df["cop_cor.IN_SAM", "delta_n"] <- "-"
  return(result_df)
}


combine_result_tables <-
  function(result_tables, statistic_names = c(
             "ME(%)", "STD(%)", "SK", "KU", "SR",
             "$\\bar{U}(10^4)$", "TO(%)", "$SR_{net}$",
             "$\\bar{U}_{net}(10^4)$", "$\\Delta_{net}$(%)"
           ),
           strategy_names = c(
             "动态-t", "动态-N", "静态-t", "静态-N", "历史",
             "固定O", "固定I"
           )) {
    #' @title 对于一个sum1/nosum1 下的表格list，整理需要的行，然后将不同gamma 表格合并
    #' @param result_tables 不同gamma 下的结果表格list。
    #' @return data.frame 一个将不同gamma 下的策略合并在一起的大的result table

    # 整理每个一个gamma 下结果表格的格式
    formatted_result_tables <- lapply(result_tables, .format_result_table,
      annualize_scale = 52
    )

    # list 中加入插入空行，然后合并所有表格
    statistic_num <- ncol(formatted_result_tables[[1]])
    tables_with_split <- formatted_result_tables
    split_row <- list(matrix(rep("", statistic_num), nrow = 1))
    for (i in (length(tables_with_split) - 1):0) {
      tables_with_split <- append(tables_with_split, split_row, i)
    }
    result_binded <- do.call("rbind", tables_with_split)

    # 更改行名、列名
    gamma_names <- paste0("$\\boldsymbol{\\gamma = ", names(result_tables), "}$")
    toal_rownames <- c()
    for (gamma_name in gamma_names) {
      toal_rownames <- c(toal_rownames, gamma_name, strategy_names)
    }
    rownames(result_binded) <- toal_rownames
    colnames(result_binded) <- statistic_names
    return(result_binded)
  }
