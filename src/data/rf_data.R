# ==================================================================
# 整理无风险收益文件的格式
# ==================================================================

suppressPackageStartupMessages({
    library(xts)
})

rf_df <- read.csv("data/raw/risk_free.csv",
    header = TRUE, sep = "\t",
    fileEncoding = "UCS-2LE", stringsAsFactors = FALSE
)

rf_xts <- xts::as.xts(rf_df[, 2:ncol(rf_df)], order.by = as.Date(rf_df[, 1]))
# 原始数据是百分比，需要转换成小数
rf_xts <- rf_xts / 100
colnames(rf_xts) <- c("Year_rf", "Day_rf", "Week_rf", "Month_rf")
saveRDS(rf_xts, "data/interim/rf_xts.Rds")
