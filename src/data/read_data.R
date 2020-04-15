#############################################################################
# 读取中间数据的接口
#############################################################################

read_fac_xts <- function(data_freq = "Week") {
  #' @title 读取data preparing 之后的xts 对象

  #' @param data_freq Default "Week" 指定哪个频率的数据，c("Week", "Month", "Day") 的其中之一
  #' @return 一个xts 对象，为保存过的因子xts 数据
  data_path <- paste0("data/interim/fac_xts_", data_freq, ".Rds")
  xts_obj <- readRDS(data_path)
  return(xts_obj)
}


read_best_arma_order <- function(data_path =
                                   "data/interim/best_arma_ssdt_Week.Rds") {
  #' @title 读取最佳arma order 数据
  #' @return 返回一个matrix, 保存每个因子的最佳arma order。
  best_arma_order <- readRDS(data_path)
  return(best_arma_order)
}
