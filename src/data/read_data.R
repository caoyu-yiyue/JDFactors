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


read_rolling_multigarchfit <-
  function(data_path = "data/interim/rolling_multigarch.Rds") {
    #' @title 读取rolling mutiplegarch fit 对象的list
    #' @param data_path 保存数据的路径
    #' @return list of uGARCHmultifit，以拟合garch 时的数据行数为names

    multigarch_list <- readRDS(data_path)
    return(multigarch_list)
  }


read_rolling_cop_rcov <-
  function(data_path = "data/interim/rolling_cop_rcov.Rds", which) {
    #' @title 读取rolling copula fit 得到的滚动一步预测方差-协方差矩阵(rcov)
    #'
    #' @param data_path str. 保存数据的路径
    #' @param which str. one of c("t_dcc", "norm_dcc", "t_static",
    #' "norm_static", "all")
    #' 指定需要返回哪种copula 模型生成的rcov。
    #' @return which == "all" 时，将返回包含四种copula rcov 的list，每个都是xts 对象
    #' which 为其他时，返回该指定的copula 生成的rcov。

    rolling_cop_rcov_list <- readRDS(data_path)
    cop_rcov <- if (which == "all") {
      rolling_cop_rcov_list
    } else {
      rolling_cop_rcov_list[[which]]
    }

    # 如果读取的数据为NULL，应该为which 参数不合法，报错。
    if (is.null(opt_weights)) {
      stop("Read NULL Data. Seems Wrong 'which' Parameter.")
    }


    return(cop_rcov)
  }


read_rolling_mean <- function(data_path = "data/interim/rolling_mean.Rds") {
  #' @title 读取rolling mean(滚动几何平均值)的 xts 对象。
  #' @param 读取数据的路径
  #' @return xts. 即保存了rolling geom mean 的xts 对象

  rolling_mean_xts <- readRDS(data_path)
  return(rolling_mean_xts)
}


read_opt_weights <- function(data_path = "data/interim/opt_weights.Rds",
                             which) {
  rolling_cop_rcov_list <- readRDS(data_path)
  opt_weights <- if (which == "all") {
    rolling_cop_rcov_list
  } else {
    rolling_cop_rcov_list[[which]]
  }

  if (is.null(opt_weights)) {
    stop("Read NULL Data. Seems Wrong 'which' Parameter.")
  }

  return(opt_weights)
}
