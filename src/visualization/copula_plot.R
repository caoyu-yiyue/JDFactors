suppressPackageStartupMessages({
  library(xts)
  library(rmgarch)
})

source("src/data/read_data.R")

# 读取cgarchfit 对象
full_cops <- read_all_cops()
cops_list <- full_cops[c("dcc_t", "dcc_norm")]


dcc_cop_cor_plot <- function(dcc_cops, colors, legend_txt = NULL, ltys = 1:6,
                             main_title = "Corrs for Copulas",
                             xlab = "Date", ylab = "Correlation") {
  #' @title 对多个dcc copulas 画rcor 序列图。每个因子对儿为一张图，不同的cop
  #' 同样的因子对儿将画在一张图上。小图们会组合为一张大图。
  #'
  #' @param dcc_cops list. 保存了dcc = TRUE 的rmgarch::cgarchfit 对象的list
  #' @param colors 图形的颜色序列，长度需与dcc_cops 一致
  #' @param legend_txt 图例的名字, not using.
  #' @param main_title 组合图片的主标题，即放在最上面的总标题
  #' @param xlab 每个图形的x 轴名称
  #' @param ylan 每个图形的y 轴名称
  #' @return NULL 但是会在默认dev 中输出图形

  # 对传入的dcc_cops list 求rcor。根据第一个rcor 找到图数、index、图名序列
  rcor_list <- lapply(dcc_cops, rcor, output = "matrix")
  n_plots <- ncol(rcor_list[[1]])
  idx <- index(rcor_list[[1]])
  plot_names <- colnames(rcor_list[[1]])

  par(mfrow = c(ceiling(n_plots / 3), 3), oma = c(0, 0, 1, 0), xpd = TRUE)
  blank_polts <- prod(par()$mfrow) - n_plots
  # 对每个图形循环（实际上是对rcor xts 的每列的idx 循环
  for (plot_i in 1:n_plots) {
    # 使用第一个rcor xts 创建plot，同时指定xlab 和ylab
    plot(
      x = idx, y = rcor_list[[1]][, plot_i], type = "l", col = colors[1],
      xlab = xlab, ylab = ylab, lty = ltys[1]
    )

    # 对第2 到后面的rcor 们，使用lines 函数直接添加到上面的plot 上
    for (mat_i in 2:length(rcor_list)) {
      lines(
        x = idx, y = rcor_list[[mat_i]][, plot_i],
        col = colors[mat_i], lty = ltys[mat_i]
      )
    }
    title(main = plot_names[plot_i], adj = 0, line = 0.8)
  }

  # 加总标题和图例
  mtext(main_title, side = 3, outer = TRUE, line = -1.5)
  if (!is.null(legend_txt)) {
    if (legend_txt == "default") {
      legend_txt <- names(dcc_cops)
    }

    # 把然后去最后加图例空图跑完
    for (i in seq_len(blank_polts)) {
      plot.new()
    }

    legend("bottomright",
      legend = legend_txt, col = colors,
      lty = ltys[seq_len(length(dcc_cops))], xpd = NA
    )
  }
  par(mfrow = c(1, 1), oma = rep.int(0, 4), xpd = FALSE)
}

dcc_cop_cor_plot(cops_list, colors = c("green", "blue"), legend_txt = "default")
