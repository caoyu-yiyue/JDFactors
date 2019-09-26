
library(rugarch)
library(copula)
library(readr)
library(dplyr)
library(xts)
library(sn)
library(foreach)

# 读取数据并转换成xts
week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))


# garch 模型的指定
var_mdl <- list(model = "gjrGARCH", garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")


###### 滚动估计garch ######
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
arch_roll_dist_df <- lapply(arch_roll, FUN = function(x) {
  as.data.frame(x, which = "density")[, 1:4]
})

# 保存roll garch 的结果
# save(arch_roll, arch_roll_dist_df, file = "data/interim/arch_roll_dist.Rda")



#### 按年滚动估计copula 系数 #####
cls <- makeCluster(4, type = "FORK")
doParallel::registerDoParallel(cls, cores = 2)
# 这里因为最后一年只使用前一年的模型，所以不必对它进行循环计算
cop_params <- foreach(i = year_end_idx[4:(length(year_end_idx) - 1)], .combine = "rbind") %dopar% {
  d_dim <- ncol(week_fac)
  d_window <- week_fac[1:i]

  pob_var <- pobs(as.matrix(d_window))
  cop_mdl <- tCopula(dim = d_dim)
  cop_fit <- fitCopula(cop_mdl, pob_var, method = "ml")

  # year_dix_last <- format(index(week_fac[year_end_idx[i]]), "%Y")
  year_idx <- index(week_fac[i + 1])
  copula_params <- xts(t(coef(cop_fit)), year_idx)

  return(copula_params)
}
stopCluster(cls)



####### merge factors dp for skewt dist #############
# 把每天的skew-t 分布的cp 参数转化成ns 包使用的dp 参数
cpDf2dpDf <- function(df) {
  # 对df 中的每一行进行cp2dp，最后返回一个xts 对象
  dp_matrix <- t(apply(df, 1, FUN = cp2dp, family = "ST"))
  return(as.xts(dp_matrix))
}
# 对arch_roll_dist_df 中的每一个data.farme 进行转换
roll_dist_dp <- lapply(arch_roll_dist_df, cpDf2dpDf)

# 把不同因子的参数merge 到一起
names(Reduce(f = function(left, right) {
  merge.xts(as.xts(left), as.xts(right))
}, roll_dist_dp))






####### 下面的先不用管 ###########

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
