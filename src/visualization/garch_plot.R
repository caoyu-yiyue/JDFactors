#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rugarch)
  library(copula)
  library(gtools)
})


plot_multi_garch <- function(garch_fits_list, which_, fac_names, main_title) {
  #' @title 对一个list 的garch fit 对象调用plot 方法，把五个图形画在一个图上
  #'
  #' @param garch_fit_list 保存多个（五个）garch fit 的list。(保存的multigarchfit@fit)
  #' @param which_ ugarchfit plot method 的which 参数对应的数字，对应画哪种图。
  #' @param fac_names 因子名，用于放在图形下标示因子的名称。顺序需与garch_fits_list
  #' 中保存的因子顺序一致
  #' @param main_title 所有图片的主标题
  #' @return NULL 不返回任何值，但是将在默认的device 输出图片。

  # 先设定多图；主标题颜色白色，从而隐藏默认标题。
  par(mfrow = c(2, 3), col.main = "white")
  # 对每个因子循环，画出对应的图，同时在左下角加上一个因子名。
  for (i in seq_len(length(fac_names))) {
    plot(garch_fits_list[[i]], which = which_)
    order <- paste0(i, ")")
    title(sub = paste(order, fac_names[i], sep = " "), adj = 0, )
  }
  # 标题变回黑色，同时在顶部加上主标题；最后恢复单图设置。
  par(col.main = "black")
  mtext(text = main_title, side = 3, outer = TRUE, line = -3)
  par(mfrow = c(1, 1))
}


.quant_depend_coef_emp <- function(quantile_x, data_pair) {
  #' @title 计算经验的quantile dependence coeffient
  #'
  #' @param quantile_x 单个数字。输入一个分位数值，用于作为threshold
  #' @param data_pair 含有两列的数据对儿，必须为已经转换为均匀分布的数据
  #' @return quantile dependece coefficent 数值

  # 对于小于0.5 的，为二者都小于quantile 的个数，除以其中一个（第二个）小于的个数
  if (quantile_x < 0.5) {
    denominator <- nrow(data_pair[data_pair[, 2] <= quantile_x])
    molecule <- nrow(data_pair[data_pair[, 1] <= quantile_x &
      data_pair[, 2] <= quantile_x])
  } else {
    # 对于大于0.5 的，则相反。
    denominator <- nrow(data_pair[data_pair[, 2] > quantile_x])
    molecule <- nrow(data_pair[data_pair[, 1] > quantile_x &
      data_pair[, 2] > quantile_x])
  }
  return(molecule / denominator)
}


.cop_implied_qdc <- function(data_pair, cop_type) {
  #' @title 计算一个copula 隐含的quantile dependence coeﬃcient 序列
  #'
  #' @param data_pair xts 对象，一个两列数据组成的数据对儿，必须是经过PIT 转换的数据
  #' @param cop_type str one of c("normal", "t")。指定copula 的类型，这里支持normal 和t
  #' @return numeric 计算过后的qdc 值。
  #' @details 这里使用传入的数据估计一个copula。
  #' 对于估计的t-copula，df 四舍五入为整数后再计算pCopula 的值。

  # 拟合一个copula
  cop <- ellipCopula(family = cop_type, dim = 2, dispstr = "un")
  cop_fit <- fitCopula(copula = cop, data = as.matrix(data_pair), method = "ml")
  # 解析参数rho 和df（df 需要约为整数），df 需要在cop_type 为t 时使用
  rho <- coef(cop_fit)[1]
  df <- if (cop_type == "t") {
    round(coef(cop_fit)[2])
  }

  # 使用估计的参数新建一个copula 对象，df 参数依然是cop_type 为t 时指定。
  cop_paramed <- ellipCopula(
    family = cop_type, dim = 2, dispstr = "un", param = rho,
    df = if (cop_type == "t") {
      df
    }
  )

  # 计算左侧的qdc 和右侧的qdc，合并到一起并输出
  left_x <- seq(from = 0.05, to = 0.5, by = 0.01)
  left_mat <- cbind(left_x, left_x)
  left_y <- pCopula(u = left_mat, copula = cop_paramed) / left_x

  right_x <- seq(from = 0.51, to = 0.95, by = 0.01)
  right_mat <- cbind(right_x, right_x)
  right_y <- (1 - 2 * right_x + pCopula(u = right_mat, copula = cop_paramed)) /
    (1 - right_x)

  qdc <- c(left_y, right_y)
  return(qdc)
}


plot_qdc_all <- function(multigarch_fit_list, facs_names) {
  #' @title 画全部数据对儿的quantile dependece coefficent 图形
  #' @param multigarch_fit_list list of ugarchfit. 保存GARCH 模型拟合结果的list。

  # 将数据转换为均匀分布
  facs_unif <- do.call(cbind, lapply(multigarch_fit_list, pit))
  colnames(facs_unif) <- facs_names
  # 将所有的因子两两配对，形成组合。这里是对col index 进行组合以保证图片顺序。
  fac_nums <- ncol(facs_unif)
  col_num_pairs <- gtools::combinations(fac_nums, 2, seq_len(fac_nums))
  # 所有图形的x 轴，一个quantile 的序列
  quantile_x_vec <- seq(from = 0.05, to = 0.95, by = 0.01)
  order <- 1

  # ======================== 对每个col pair 画图的函数 ======================= #
  .plot_for_one_data_pair <- function(col_num_pair) {
    #' @title 该函数将被用于每个col 对儿，在这个函数中进行画图工作
    #' @param col_num_pair 输入col 序号的pair
    #' @return NULL 不返回任何值，但在标准dev 中进行画图

    data_pair <- facs_unif[, col_num_pair]
    # 对每个quantile（即X轴）应用计算quantile dependence coef 的函数，得到Y轴
    coef_y_vec_emp <- sapply(quantile_x_vec,
      .quant_depend_coef_emp,
      data_pair = data_pair
    )
    norm_cop_qdc <- .cop_implied_qdc(data_pair = data_pair, cop_type = "normal")
    t_cop_qdc <- .cop_implied_qdc(data_pair = data_pair, cop_type = "t")
    # #### 可以在这里继续计算更多的y，然后画到图里 #### #

    # 开始画图
    matplot(
      x = quantile_x_vec, y = cbind(coef_y_vec_emp, norm_cop_qdc, t_cop_qdc),
      type = c("l", "l", "l"),
      lty = c(1, 4, 5),
      xlab = "Quantile", ylab = "Quantile Dependence"
    )
    # 画0.5 处的竖线
    abline(v = 0.5, lty = "dotted")

    # 所有的图画完之后，加入下标以显示因子名
    name_pair <- paste(colnames(data_pair), collapse = " - ")
    order_txt <- paste0("(", order, ")")
    subscript <- paste(order_txt, name_pair)
    title(sub = subscript, adj = 0)
    order <<- order + 1
  }
  # ========================================================================= #

  # 对因子列的index 配对应用以上画图函数，画出图形
  # 首先规定layout，其中最下方为legend 留一整行（四个13）
  layout_mat <- matrix(c(1:12, 13, 13, 13, 13),
    nrow = 4, ncol = 4, byrow = TRUE
  )
  layout(mat = layout_mat, heights = c(rep(0.3, 3), 0.1))
  apply(col_num_pairs, 1, .plot_for_one_data_pair)

  # 跳过空图
  par(mar = c(0.5, 0, 0, 0)) # 这里是规定图片的margin，防止legend 出界
  for (i in 1:3) {
    plot.new()
  }
  legend("center",
    inset = 0,
    legend = c("Data", "norm-Copula", "t-Copula"),
    lty = c(1, 4, 5), col = 1:3, xpd = NA, horiz = TRUE,
  )
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}
