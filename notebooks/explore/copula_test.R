library(rugarch)
library(copula)
library(readr)
library(dplyr)
library(xts)
library(sn)

# 读取数据并转换成xts
week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))


# garch 模型的指定
var_mdl <- list(model = "gjrGARCH", garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")


# 滚动估计garch
year_end_idx <- endpoints(week_fac, on = "year")
cluster <- makeCluster(4, type = "FORK")
rolling_arch <- function(fac_ts) {
  return(ugarchroll(archSpec, fac_ts,
    n.ahead = 1, n.start = year_end_idx[4],
    refit.every = 1, refit.window = "recursive",
    calculate.VaR = FALSE, keep.coef = FALSE, cluster = cluster
  ))
}
arch_roll <- lapply(week_fac, rolling_arch)
stopCluster(cluster)
arch_roll_dist_df <- lapply(arch_roll, FUN = as.data.frame, which = "density")

# 保存roll garch 的结果
# save(arch_roll, arch_roll_dist_df, file = "data/interim/arch_roll_dist.Rda")


# 使用前三年的数据计算的t_copula，然后使用预测的第一天计算的因子分布得出多维联合分布
test_windwo_1 <- week_fac[1:year_end_idx[4]]
cor(test_windwo_1, method = "kendall")

pob_var <- pobs(as.matrix(test_windwo_1))
cop_mdl <- tCopula(dim = 3)
cop_fit <- fitCopula(cop_mdl, pob_var, method = "ml")

rho <- coef(cop_fit)[1]
df <- coef(cop_fit)[2]

margin_params <- sapply(arch_roll_dist_df, function(x) {
  cp2dp(as.numeric(x[1, 1:4]), family = "ST")
})

mv_dist <- mvdc(
  copula = tCopula(rho, dim = 3, df = df), margins = c("st", "st", "st"),
  paramMargins = list(
    dp = margin_params[, "premium"], dp = margin_params[, "smb"],
    dp = margin_params[, "hml"]
  )
)
rMvdc(100, mv_dist)
