#!/usr/bin/env Rscript

#############################################################################
# 计算facs 的52 周历史滚动几何平均收益率的功能
# 单独运行脚本时，需要输入一个参数：保存数据的路径
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(zoo)
  library(EnvStats)
  library(optparse)
})

source("src/data/read_data.R")
source("src/config.R")


.geom_mean_return <- function(x, na_rm = FALSE) {
  #' @title 输入一组收益率，计算得到它们的几何平均收益率
  #' @param x numeric. 一组收益率vector
  #' @param na_rm 是否移除NA 值
  #' @return scalar. 几何收益率

  # 使用EnvStats::geoMean，以输出它定义的错误信息。
  return(EnvStats::geoMean(1 + x, na_rm) - 1)
}


rolling_mean_main <- function() {
  #' @title 计算因子的滚动历史几何收益率的主函数。生成xts 保存到命令行传入的路径中。
  #' @details 注意：每个日期计算它之前52 期（而不包含当期）的几何收益率。
  #' @return NULL 但将会保存一个xts 对象到命令行传入的第一个参数说明的路径中。

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

  # 读取数据，并且只使用rolling copula 时期的数据，避免多余计算
  facs_xts <- read_fac_xts(data_freq = data_freq)
  in_sapmle_end <- in_sample_yearend_row(facs_xts, IN_SAMPLE_YEARS[data_freq])
  window_len <- 52
  # 须转换为zoo 对象以使用rollapply。同时所需数据为in_sample_end - window_len 后 + 1
  data_for_use <- as.zoo(
    facs_xts[(in_sapmle_end - window_len + 1):nrow(facs_xts), ]
  )

  # 滚动计算几何平均收益率，每列单独计算
  rolling_geom_mean <- rollapply(
    data = data_for_use, width = list(seq(-window_len, -1)),
    FUN = .geom_mean_return, by.column = TRUE
  )
  rolling_geom_mean_xts <- xts::as.xts(rolling_geom_mean)

  # 保存到输入的路径当中。
  saveRDS(rolling_geom_mean_xts, file = opts[["out_put"]])
}


if (sys.nframe() == 0) {
    rolling_mean_main()
}
