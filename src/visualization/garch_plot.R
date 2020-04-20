#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rugarch)
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
