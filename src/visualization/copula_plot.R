################################################################################
# 针对rmgarch::cGARCHfit 对象的部分绘图函数
################################################################################

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
  cops_num <- length(dcc_cops)
  rcor_list <- lapply(dcc_cops, rcor, output = "matrix")
  n_plots <- ncol(rcor_list[[1]])
  idx <- index(rcor_list[[1]])
  plot_names <- colnames(rcor_list[[1]])

  layout_mat <- if (!is.null(legend_txt)) {
    matrix(c(1:12, 13, 13, 13, 13), nrow = 4, ncol = 4, byrow = TRUE
  ) } else {
    matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
  }

  layout(mat = layout_mat, heights = c(rep(0.3, 3), 0.1))
  blank_polts <- length(unique(as.numeric(layout_mat))) - n_plots
  # 对每个图形循环（实际上是对rcor xts 的每列的idx 循环
  for (plot_i in 1:n_plots) {
    # 使用第一个rcor xts 创建plot，同时指定xlab 和ylab
    plot(
      x = idx, y = rcor_list[[1]][, plot_i], type = "l", col = colors[1],
      xlab = xlab, ylab = ylab, lty = ltys[1]
    )

    # 对第2 到后面的rcor 们，使用lines 函数直接添加到上面的plot 上
    if (cops_num > 1) {
      for (mat_i in 2:length(rcor_list)) {
        lines(
          x = idx, y = rcor_list[[mat_i]][, plot_i],
          col = colors[mat_i], lty = ltys[mat_i]
        )
      }
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
    par(mar = c(0.5, 0, 0, 0)) # 这里是规定图片的margin，防止legend 出界
    for (i in seq_len(blank_polts)) {
      plot.new()
    }

    legend("center",
      legend = legend_txt, col = colors,
      lty = ltys[seq_len(length(dcc_cops))], xpd = NA,
      horiz = TRUE
    )
  }
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}
