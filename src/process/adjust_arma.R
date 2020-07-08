#!/usr/bin/env Rscript

#############################################################################
# 手动调整best arma 数据框中的值，使得数据通过（自相关）检验
#############################################################################

source("src/data/read_data.R")
best_arma <- read_best_arma_order(which = "origin")
best_arma[, "SMB"] <- c(3, 3)
saveRDS(best_arma, file = "data/interim/best_arma_adjusted.Rds")
