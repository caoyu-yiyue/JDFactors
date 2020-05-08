#############################################################################
# 读取中间数据的接口
#############################################################################

read_fac_xts <- function(data_freq = "Week") {
  #' @title 读取data preparing 之后的xts 对象

  #' @param data_freq Default "Week" 指定哪个频率的数据，c("Week", "Month", "Day") 的其中之一
  #' @return 一个xts 对象，为保存过的因子xts 数据

  all_facs_xts_list <- readRDS("data/interim/facs_xts.Rds")
  xts_obj <- all_facs_xts_list[[data_freq]]
  return(xts_obj)
}


read_best_arma_order <- function(data_path =
                                   "data/interim/best_arma_ssdt_Week.Rds") {
  #' @title 读取最佳arma order 数据
  #' @return 返回一个matrix, 保存每个因子的最佳arma order。
  best_arma_order <- readRDS(data_path)
  return(best_arma_order)
}


read_multi_garch_fit <- function(data_path =
                                   "data/interim/multi_garch_mdl.Rds") {
  #' @title 读取多变量multi garch fit 对象
  #' @return rugarch::uGARCHmultifit 即数据没列进行fit 之后的multifit 对象
  multi_garch_fit_obj <- readRDS(data_path)
  return(multi_garch_fit_obj)
}


read_all_cops <- function(data_path = "data/processed/all_cops.Rds") {
  #' @title 读取针对全时常数据计算的所有copula 模型的list
  # ‘ @return 保存所有copula 的list
  all_cops <- readRDS(data_path)
  return(all_cops)
}


in_sample_yearend_row <- function(data, in_sample_year) {
  #' @title 针对数据的in_sample_year 的年数，找到in_sample 的最后一行行号
  #' @param data xts 对象，即数据本身
  #' @param in_sample_year int，指定in sample year 的年数
  #' @return int. in sample 数据的最后一行行号
  year_endponits <- xts::endpoints(data, on = "year")[-1]
  in_sample_end_row <- year_endponits[in_sample_year]
  return(in_sample_end_row)
}
