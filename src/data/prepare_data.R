# ==================================================================================
# 从数据库中下载的原始数据开始，读取原始数据、形成所需的数据格式等准备工作
# ==================================================================================
library(xts)
library(ISOweek)

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
  add_weekday <- paste0(substr(week_col, 1, 5), "W", substr(week_col, 6, 7), "-5")
  date_col <- ISOweek2date(add_weekday)

  df_with_week[[week_col_name]] <- date_col
  return(df_with_week)
}


if (!interactive()) {
  week_raw <- read_raw(data_freq = "Week")
  week_to_date <- parse_year_week(week_raw)
}
