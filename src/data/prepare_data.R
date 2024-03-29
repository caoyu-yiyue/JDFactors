# ==============================================================================
# 从数据库中下载的原始数据开始，读取原始数据、形成所需的数据格式等准备工作
# ==============================================================================
suppressPackageStartupMessages({
  library(ISOweek)
  library(xts)
})
read_raw <- function(data_freq) {
  # 读取从数据库下载到的原数据。
  # Arguments: data_freq 哪种数据频率的原数据, c("Week", "Month", "Day") 的其中之一
  # Value: data.frame 原数据的data.frame 对象

  # 使用传入的data_freq 组成数据路径
  path <- paste0("data/raw/Fivefac", data_freq, ".csv")
  # 读取数据
  raw_d <- read.csv(path,
    header = TRUE, sep = "\t",
    fileEncoding = "UCS-2LE", stringsAsFactors = FALSE
  )
  return(raw_d)
}


parse_year_week <- function(df_with_week, week_col_name = "TradingWeek") {
  # 格式化year-week 的日期列为当周周五的date
  # Arguments:
  #   df_with_week: 带有week 列的data.frame
  #   week_col_name: week 列的列名
  # Value:
  #   data.frame 将week 列转换为date 的data.frame

  week_col <- df_with_week[[week_col_name]]
  # 将week 列按照IOSweek2date 所需的形式，转换为"%Y-W%V-%u", 这里直接设置为星期五
  add_weekday <- paste0(
    substr(week_col, 1, 5), "W",
    substr(week_col, 6, 7), "-5"
  )
  date_col <- ISOweek2date(add_weekday)

  df_with_week[[week_col_name]] <- date_col
  return(df_with_week)
}


fac_df_to_xts <- function(fac_df) {
  # 将因子的data.frame 对象转换为xts 对象
  # Arguments:
  #   fac_df: 保存有因子（第一列为时间列，类型为Date）的data.frame
  # Values:
  #   xts 对象，保存有原数据的xts 对象。

  # 首先将列名中的1 和Trading 去掉；将RiskPremium 改为Market
  colnames(fac_df) <- gsub("1|Trading", "", colnames(fac_df))
  colnames(fac_df) <- sub("RiskPremium", "Market", colnames(fac_df))
  # 转换xts 对象并返回
  xts_obj <- as.xts(fac_df[, -1], order.by = fac_df[, 1])
  return(xts_obj)
}


prepare_data_main <- function() {
  # 进行数据预备的主函数，把Week, Day, Month 的原始csv 数据都保存为xts 对象。
  # 三个xts 对象将组成一个list，name 分别为Week, Day, Month.
  # 保存路径名为data/interim/facs_xts.Rds

  all_xts_obj <- list()
  for (freq in c("Week", "Month", "Day")) {
    raw_data <- read_raw(data_freq = freq)

    if (freq == "Week") {
      raw_data <- parse_year_week(raw_data)
    } else if (freq == "Month") {
      # Month 数据将TradingMonth 列转为yearmon 类型
      raw_data[["TradingMonth"]] <- zoo::as.yearmon(raw_data[["TradingMonth"]])
    } else if (freq == "Day") {
      # Day 数据将TradingDate 列转为Date 类型
      raw_data[["TradingDate"]] <- as.Date(raw_data[["TradingDate"]])
    }

    # 转换为xts 对象并加入到总list 中
    xts_obj <- fac_df_to_xts(raw_data)
    all_xts_obj[[freq]] <- xts_obj
  }

  saveRDS(all_xts_obj, file = "data/interim/facs_xts.Rds")
}


if (sys.nframe() == 0) {
  prepare_data_main()
}
