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


read_rf_xts <- function(data_freq = "Week") {
  #' @title 读取无风险收益率
  #' @param data_freq 数据频度，one of c("Week", "Month", "Day", "Year", "all")
  #' all 将返回所有的数据
  #' @return xts 对象。data_freq == "all" 将返回四列；else 返回一列

  all_rf_xts <- readRDS("data/interim/rf_xts.Rds")
  rf_xts <- if (data_freq == "all") {
    all_rf_xts
  } else {
    all_rf_xts[, paste0(data_freq, "_", "rf")]
  }

  return(rf_xts)
}


read_best_arma_order <- function(which) {
  #' @title 读取最佳arma order 数据
  #' @return 返回一个matrix, 保存每个因子的最佳arma order。
  #'
  if (which == "origin") {
    data_path <- "data/interim/best_arma_ssdt_Week.Rds"
  } else if (which == "adjusted") {
    data_path <- "data/interim/best_arma_adjusted.Rds"
  } else {
    stop("Wrong parameter 'which', must one of 'origin' and 'adjusted'.")
  }
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

################################################################################
# 与rolling fit 相关的数据读取 部分
################################################################################

in_sample_yearend_row <- function(data, in_sample_year) {
  #' @title 针对数据的in_sample_year 的年数，找到in_sample 的最后一行行号
  #' @param data xts 对象，即数据本身
  #' @param in_sample_year int，指定in sample year 的年数
  #' @return int. in sample 数据的最后一行行号
  year_endponits <- xts::endpoints(data, on = "year")[-1]
  in_sample_end_row <- year_endponits[in_sample_year]
  return(in_sample_end_row)
}


read_rolling_best_arma <-
  function(base_name = "data/interim/rolling_best_arma",
           data_freq = "Week", extension = ".Rds",
           which = "best_arma") {
    #' @title 读取rolling best arma 的数据
    #' @param data_freq str. 数据频度，默认"Week"
    #' @param which str. one of c("best_arma", "all_bics") 指定读取最佳arma
    #' 还是计算的所有的bics。

    data_path <- paste0(base_name, "_", data_freq, extension)
    best_armas <- readRDS(data_path)
    if (which == "best_arma") {
      return(best_armas[["best_armas"]])
    } else if (which == "all_bics") {
      return(best_armas[["all_bics"]])
    } else {
      stop("Wrong 'which' Parameter. Must one of 'best_arma' and 'all_bics'")
    }
  }


read_rolling_multigarchfit <-
  function(base_name = "data/interim/rolling_multigarch", data_freq = "Week",
           extension = ".Rds") {
    #' @title 读取rolling mutiplegarch fit 对象的list
    #' @param data_path 保存数据的路径
    #' @return list of uGARCHmultifit，以拟合garch 时的数据行数为names

    data_path <- paste0(base_name, "_", data_freq, extension)
    multigarch_list <- readRDS(data_path)
    return(multigarch_list)
  }


read_rolling_cop_rcov <-
  function(which, base_name = "data/interim/rolling_cop_rcov",
           data_freq = "Week", extension = ".Rds") {
    #' @title 读取rolling copula fit 得到的滚动一步预测方差-协方差矩阵(rcov)
    #'
    #' @param data_path str. 保存数据的路径
    #' @param which str. one of c("t_dcc", "norm_dcc", "t_static",
    #' "norm_static", "all")
    #' 指定需要返回哪种copula 模型生成的rcov。
    #' @return which == "all" 时，将返回包含四种copula rcov 的list，每个都是xts 对象
    #' which 为其他时，返回该指定的copula 生成的rcov。

    data_path <- paste0(base_name, "_", data_freq, extension)
    rolling_cop_rcov_list <- readRDS(data_path)
    cop_rcov <- if (which == "all") {
      rolling_cop_rcov_list
    } else {
      rolling_cop_rcov_list[[which]]
    }

    # 如果读取的数据为NULL，应该为which 参数不合法，报错。
    if (is.null(cop_rcov)) {
      stop("Read NULL Data. Seems Wrong 'which' Parameter.")
    }


    return(cop_rcov)
  }


read_fixed_cor_rcov <-
  function(which, base_name = "data/interim/fixed_cor2cov",
           data_freq = "Week", extension = ".Rds") {
    #' @title 读取使用样本内（外）cor matrix 得到的滚动一步预测方差-协方差矩阵(rcov)
    #'
    #' @param data_path str. 保存数据的路径
    #' @param which str. one of c("fixed_cor.IN_SAM", "fixed_cor.OUT_SAM"")
    #' 指定需要返回哪种copula 模型生成的rcov。
    #' @return which == "all" 时，将返回包含样本内（外）cor 计算的两个rcov 的list，每个都是xts 对象
    #' which 为其他时，返回该指定的cor matrix 生成的rcov。

    data_path <- paste0(base_name, "_", data_freq, extension)
    fixed_cor_gen_rcov_list <- readRDS(data_path)
    fixed_cor_rcov <- if (which == "all") {
      fixed_cor_gen_rcov_list
    } else {
      fixed_cor_gen_rcov_list[[which]]
    }

    # 如果读取的数据为NULL，应该为which 参数不合法，报错。
    if (is.null(fixed_cor_rcov)) {
      stop("Read NULL Data. Seems Wrong 'which' Parameter.")
    }

    return(fixed_cor_rcov)
  }


read_rolling_mean <- function(base_name = "data/interim/rolling_mean",
                              data_freq = "Week", extension = ".Rds") {
  #' @title 读取rolling mean(滚动几何平均值)的 xts 对象。
  #' @param 读取数据的路径
  #' @return xts. 即保存了rolling geom mean 的xts 对象

  data_path <- paste0(base_name, "_", data_freq, extension)
  rolling_mean_xts <- readRDS(data_path)
  return(rolling_mean_xts)
}


read_opt_weights <- function(which, base_name = "data/interim/opt_weights",
                             data_freq = "Week", extension = ".Rds") {
  #' @title 读取滚动最优化求得的权重值
  #' @param data_path str. 保存数据的路径
  #' @param which three item charactor vector.
  #' 第一个元素指定是否需要sum1，即one of c("sum1", "no_sum1")
  #' 第二个元素指定生成数据时的gamma，one of c(3, 8, 20, 50)；
  #' 第三个元素为生成数据时的copula type,
  #' one of c("t_dcc", "norm_dcc", "t_static", "norm_static", "all")
  #' 指定需要返回哪种copula 模型计算的最优权重。
  #' 或者可以传入一个"all"，返回包含所有结果的整个list。
  #' 另外，第三或二三个项目可以为空，将返回第二级或第一级的list。
  #' @return which == "all" 时，将返回包含四种copula 最优权重的list，每个都是xts 对象
  #' which 为其他时，返回该指定的copula 计算所得的权重。


  data_path <- paste0(base_name, "_", data_freq, extension)
  rolling_cop_rcov_list <- readRDS(data_path)
  opt_weights <- if (which[1] == "all") {
    rolling_cop_rcov_list
  } else {
    result <- rolling_cop_rcov_list
    for (i in which) {
      result <- result[[i]]
    }
    result
  }

  if (is.null(opt_weights)) {
    stop("Read NULL Data. Seems Wrong 'which' Parameter Pair.")
  }

  return(opt_weights)
}


read_port_ret <- function(base_name = "data/processed/port_ret",
                          data_freq = "Week", extension = ".Rds") {
  #' @title 读取组合收益率list

  data_path <- paste0(base_name, "_", data_freq, extension)
  return(readRDS(data_path))
}


read_result_tables <-
  function(which, base_name = "data/processed/result_tables",
           data_freq = "Week", extension = ".Rds") {
    #' @title 读取脚本计算的滚动结果
    #' @param which str. one of c("all", "sum1", "no_sum1") 读取哪种结果
    #' @param data_freq str. 指定数据频度
    #' @return list of list if which == "all", list else.
    #' 最底层的项目为一个gamma 下的结果list。

    data_path <- paste0(base_name, "_", data_freq, extension)
    result <- if (which == "all") {
      readRDS(data_path)
    } else {
      readRDS(data_path)[[which]]
    }

    if (is.null(result)) {
      stop("'Get NULL data, Seems Wrong 'which' Parameter.")
    }
    return(result)
  }
