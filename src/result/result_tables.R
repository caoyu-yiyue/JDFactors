#!/usr/bin/env Rscript

#############################################################################
# 输出投资结果的表格时，计算所需指标的函数
#############################################################################

suppressPackageStartupMessages({
  library(xts)
  library(PerformanceAnalytics)
})

source("src/data/read_data.R")


basic_statistics_table <-
  function(ret_xts, data_freq,
           statistics_name =
             c("Return", "Std", "Skew", "Kurt", "Sharpe Ratio")) {
    #' @title 计算均值、方差、偏度、峰度、Sharp Ratio 等基础统计量的table
    #'
    #' @param ret_xts xts 对象。组合收益率，每列为一种策略，每行为一日
    #' @param data_freq chatactor, one of c("Week", "Day", "Month")，即收益率的频率
    #' @param statistics_name 单独指定统计量的名字。
    #' @return data.frame 每列为策略名，每行为统计量的一个data.frame

    scale <- switch(data_freq, "Week" = 52, "Day" = 252, "Month" = 12)
    annual_table <- table.AnnualizedReturns(ret_xts, scale = scale)
    skew <- skewness(ret_xts)
    kurt <- kurtosis(ret_xts)

    result_table <- rbind(annual_table[1:2, ], skew, kurt, annual_table[3, ])
    rownames(result_table) <- statistics_name
    return(result_table)
  }


utility_value <- function(portfolio_ret, risk_coef) {
  #' @title 计算某一种策略下组合收益的效用函数值
  #'
  #' @param portfolio_ret xts 对象，组合收益率。
  #' @param risk_coef numeric scalar. 风险厌恶系数。
  #' @return numeric scalar. 效用函数值

  port_mean <- mean(portfolio_ret)
  variance <- var(portfolio_ret)
  u_value <- port_mean - risk_coef / 2 * variance
  return(u_value)
}


single_strategy_turnover <- function(strategy_name, port_ret_xts,
                                     weights_list, facs_xts, rf_xts = NULL) {
  #' @title 计算某一个策略的平均turnover（换手率）
  #'
  #' @param strategy_name charactor, 策略名。one of c("t_dcc", "norm_dcc",
  #' "t_static", "norm_static", "fixed_cor.IN_SAM", "fixed_cor.OUT_SAM",
  #' "static_benchmark")
  #' @param port_ret_xts xts 对象，每一列是一种策略的组合收益率
  #' @param weights list，每个item 是某种策略的权重，names 为如上的strategy_name 的集合
  #' @param facs_xts xts 对象，因子每日收益率xts。
  #' @details weights 传入Begin of Period Weights，即每日的起始权重与当日的因子收益相乘得出当日的组合收益率。
  #' @return numeric scalar. 日期平均的换手率。

  # 得到因子数量，解析出策略对应的组合收益和权重值
  fac_num <- ncol(facs_xts)
  port_ret <- port_ret_xts[, strategy_name]
  weights_xts <- weights_list[[strategy_name]]

  # 如果rf_xts 不是NULL，那么将其合并到facs_xts 中
  facs_xts <- merge.xts(facs_xts, rf_xts, join = "left")

  # 每一列都除以1 加 port_ret
  end_weights <- sweep(weights_xts * (1 + facs_xts), 1,
    (1 + port_ret),
    FUN = "/"
  )
  # 每天的起始weight 减去前一天的结束weights
  weights_diff <- weights_xts - lag(end_weights)

  # 每日相加得到当日的换手率，取平均得到整体的平均换手率
  turnovers <- rowSums(abs(weights_diff[, 1:fac_num]))
  turnover_mean <- mean(turnovers, na.rm = TRUE)
  return(turnover_mean)
}


result_table_main <-
  function(port_ret_list, weights_list, facs_xts, rf_xts = NULL,
           data_freq = "Week",
           statistic_names = c("ME", "VOL", "SK", "KU", "SR", "U", "TO"),
           strategy_names = NULL) {
    #' @title 计算组合收益率统计表的主函数，将不同的计算结合在一起，再将不同的gamma 值的表合并为list
    #'
    #' @param port_ret_list portfolio return list. list names 为不同的风险偏好系数gamma，
    #' 每个项目为策略收益的xts 对象，每列为策略名，每行为时间。
    #' @param weights_list list. 为二级list，一级名gamma，二级名为策略名。每个项目为因子权重xts，
    #' 每列为因子名，每行为时间。
    #' @param facs_xts xts 对象。因子收益的xts 对象。
    #' @param rf_xts xts 对象。无风险收益xts 对象。
    #' @param data_freq charactor. 因子数据的频率。默认为"Week"，用于计算年化利率时确定scale
    #' @param statistic_names charactor. 用于单独指定统计量的名字。
    #' @param strategy_names 策略名字，默认为NULL。如果是NULL 时将使用之前代码内部的名字。
    #' @details port_ret_list 和weights_list 都接受一级名为gamma 的list，
    #' 前者每个项目为xts，列名指定策略名；后者每个项目为list，names 为策略名。
    #' 两个list 其实都是之间保存的结果文件中的sum1 / no_sum1 下的子list，
    #' 所以读取数据并传入该函数时，请自行解析出sum1 / no_sum1 下的独立项目并传入。
    #' 这样方便简化迭代深度，同时比较方便debug 两种情况中不同的权重列数可能带来的问题。
    #' @return list of matrix. names 为gamma，每个项目为该gamma 下的不同策略的统计量matrix，行名为策略名，
    #' 列名为统计量名。

    # 获取到所有的风险偏好系数，准备一个保存所有结果的空list
    risk_coefs <- names(port_ret_list)
    tables_list <- list()

    # 对每个risk_coef 循环，不使用apply 因为需要获取到risk_coef 的值
    for (risk_coef in risk_coefs) {
      # 使用risk_coef 解析得到对应的组合收益率xts，对应的权重list，策略名序列
      port_ret_xts <- port_ret_list[[risk_coef]]
      single_gamma_weights_list <- weights_list[[risk_coef]]
      strategies_vec <- names(single_gamma_weights_list)

      # 计算基础的统计量
      basic_statistics <- basic_statistics_table(
        ret_xts = port_ret_xts,
        data_freq = data_freq
      )
      # 对收益率xts 每列计算效用值
      utilities <- apply(port_ret_xts, 2,
        FUN = utility_value,
        risk_coef = as.numeric(risk_coef)
      )
      # 对每个策略计算名，计算平均换手率
      turnovers <- sapply(strategies_vec,
        FUN = single_strategy_turnover,
        port_ret_xts = port_ret_xts, weights_list = single_gamma_weights_list,
        facs_xts = facs_xts, rf_xts = rf_xts
      )

      # 合并几种统计结果，根据传入参数进行结果的重命名
      single_gamma_result <- t(rbind(basic_statistics, utilities, turnovers))
      if (!is.null(strategy_names)) {
        rownames(single_gamma_result) <- strategy_names
      }
      colnames(single_gamma_result) <- statistic_names

      # 追加到结果list 中
      tables_list[[risk_coef]] <- single_gamma_result
    }

    return(tables_list)
  }
