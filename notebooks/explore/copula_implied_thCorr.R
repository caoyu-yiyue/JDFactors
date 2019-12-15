library(rlist)
cop_implied_thCorr <- function(two_cols_df, cop_mdl) {
  # 传入一个两列的df 和一个copula model


  bottom_left <- data.frame(matrix(nrow = 0, ncol = 2))
  top_right <- data.frame(matrix(nrow = 0, ncol = 2))

  # cop_mdl <- ellipCopula(family = "norm", dim = 2)

  for (u in seq(0, 0.5, by = 0.01)) {
    qua_col_1 <- quantile(two_cols_df[, 1], u)
    qua_col_2 <- quantile(two_cols_df[, 2], u)
    extremes_xts <- two_cols_df[two_cols_df[, 1] < qua_col_1 & two_cols_df[, 2] < qua_col_2]

    # 如果行小于20 则下一循环
    if (nrow(extremes_xts) < 20) {
      next()
    }

    pob_extremes <- pobs(as.matrix(extremes_xts))
    # cop_mdl <- ellipCopula(family = "norm", dim = 2)
    fit <- fitCopula(copula = cop_mdl, pob_extremes, method = "itau")
    rho <- coef(fit)[1]

    bottom_left <- rbind(bottom_left, c(u = u, rho = rho))
  }
  colnames(bottom_left) <- c("u", "rho")

  for (u in seq(0.51, 1, by = 0.01)) {
    qua_col_1 <- quantile(two_cols_df[, 1], u)
    qua_col_2 <- quantile(two_cols_df[, 2], u)
    extremes_xts <- two_cols_df[two_cols_df[, 1] > qua_col_1 & two_cols_df[, 2] > qua_col_2]

    # 如果行小于20 则下一循环
    if (nrow(extremes_xts) < 20) {
      break()
    }

    pob_extremes <- pobs(as.matrix(extremes_xts))
    # cop_mdl <- ellipCopula(family = "norm", dim = 2)
    fit <- fitCopula(copula = cop_mdl, pob_extremes, method = "itau")
    rho <- coef(fit)[1]

    top_right <- rbind(top_right, c(u = u, rho = rho))
  }
  colnames(top_right) <- c("u", "rho")

  return(list(bottom_left, top_right))
}


fac_names <- colnames(week_fac)
col_pairs <- combn(fac_names, 2)

norm_cop_mdl <- ellipCopula(family = "norm", dim = 2)

norm_thCorr_list <- list()
for (col_num in 1:ncol(col_pairs)) {
  single_col_pair <- col_pairs[, col_num]
  threshold_corr <- cop_implied_thCorr(week_fac[, single_col_pair], cop_mdl = norm_cop_mdl)

  norm_thCorr_list <- list.append(norm_thCorr_list, threshold_corr)
}

pair_names <- apply(col_pairs, MARGIN = 2, paste, collapse = "&")
names(norm_thCorr_list) <- pair_names

















