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
  if (quantile_x <= 0.5) {
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


plot_qdc_all <- function(facs_xts) {
  #' @title 画全部数据对儿的quantile dependece coefficent 图形
  #' @param facs_xts xts 对象，即保存所有因子数据的xts对象

  # 将数据转换为均匀分布
  facs_unif <- copula::pobs(facs_xts)
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
    # #### 可以在这里继续计算更多的y，然后画到图里 #### #

    # 开始画图
    plot(quantile_x_vec, coef_y_vec_emp,
      type = "l",
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
  par(mfrow = c(3, 4))
  apply(col_num_pairs, 1, .plot_for_one_data_pair)
  par(mfrow = c(1, 1))
}
