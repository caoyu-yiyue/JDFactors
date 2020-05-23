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
colnames(rf_xts) <- c("rf_data", "day", "week", "month")
saveRDS(rf_xts, "data/interim/rf_xts.Rds")
