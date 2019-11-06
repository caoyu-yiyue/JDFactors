
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
year_end_idx <- endpoints(week_fac, on = "year")


# garch 模型的指定
var_mdl <- list(model = "gjrGARCH", garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")


###### 滚动估计garch ######
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
  as.data.frame(x, which = "density")
})

# 保存roll garch 的结果
# save(arch_roll, arch_roll_dist_df, file = "data/interim/arch_roll_dist.Rda")
load(file = "data/interim/arch_roll_dist.Rda")

arch_roll_dist_df <- lapply(arch_roll_dist_df, function(x) {
  x[, 1:4]
})


#### 按年滚动估计copula 系数 #####
cls <- makeCluster(4, type = "FORK")
doParallel::registerDoParallel(cls, cores = 2)
# 这里因为最后一年只使用前一年的模型，所以不必对它进行循环计算
cop_params <- foreach(i = year_end_idx[4:(length(year_end_idx) - 1)], .combine = "rbind") %dopar% {
  d_dim <- ncol(week_fac)
  d_window <- week_fac[1:i]

  pob_var <- pobs(as.matrix(d_window))
  cop_mdl <- tCopula(dim = d_dim, dispstr = "un")
  cop_fit <- fitCopula(cop_mdl, pob_var, method = "ml")

  # year_dix_last <- format(index(week_fac[year_end_idx[i]]), "%Y")
  year_idx <- index(week_fac[i + 1])
  copula_params <- xts(t(coef(cop_fit)), order.by = year_idx)

  return(copula_params)
}
stopCluster(cls)



####### merge factors dp for skewt dist #############
# 把每天的skew-t 分布的cp 参数转化成ns 包使用的dp 参数
cpDf2dpDf <- function(df) {
  # 对df 中的每一行进行cp2dp，最后返回一个xts 对象
  dp_matrix <- t(apply(df, 1, FUN = cp2dp, family = "ST"))

  # 确保index 为Date 对象，不包括时间
  idx <- as.Date(rownames(dp_matrix))
  return(xts(dp_matrix, order.by = idx))
}
# 对arch_roll_dist_df 中的每一个data.farme 进行转换
roll_dist_dp <- lapply(arch_roll_dist_df, cpDf2dpDf)


# 首先对每个因子循环，将其分布参数的名字前加上其自己因子名的前缀
merge_margin_dist_params <- function(margin_params_list) {
  # margin_params_list: 保存边际分布参数的list 对象，其中每一个分量是生成一个 margin dist 的参数
  # 返回一个xts 对象。合并了list 中每一个margin 的dist parameter，以其在list 中的name 为前缀。
  fac_names <- names(margin_params_list) # 一个名字的charactor vector
  for (fac in fac_names) {
    dist_param_names <- names(margin_params_list[[fac]])
    dist_param_names_with_fac <- paste(fac, dist_param_names, sep = "_")
    names(margin_params_list[[fac]]) <- dist_param_names_with_fac
  }
  # 把不同因子的分布参数merge 到一起
  merged_margin_dist_params <- Reduce(f = function(left, right) {
    merge.xts(as.xts(left), as.xts(right))
  }, x = margin_params_list)

  return(merged_margin_dist_params)
}
margin_params_xts <- merge_margin_dist_params(margin_params_list = roll_dist_dp)


####### merge magin parameters and copula parameters ###########
merged_params <- merge(margin_params_xts, cop_params, join = "outer")
# 每一年使用初年估计的copula 参数，所以使用其向前填充
merged_params <- na.locf(merged_params)
save(merged_params, file = "data/interim/merged_params.Rda")


######### 计算最佳权重 ###############
FAC_NAMES <- names(arch_roll_dist_df)
randomNum_tCop_stMargin <- function(current_row, seed = 101, n = 100000) {
  # 根据保存边际分布参数与copula 参数的行，生成联合分布的随机数
  # current_row: 保存参数的行。这里为一个单行的xts 对象
  # seed: 生成随机数的seed，默认为101。
  # n: 需要的随机数的个数，默认为100000。
  
  # 每一行单独进行操作，所以不需要index。否则后续会把index 认为是多的一列。
  current_row <- coredata(current_row)

  # 指定t_copula 对象
  t_cop_row <- tCopula(
    param = current_row[, c("rho.1", "rho.2", "rho.3")],
    dim = length(FAC_NAMES), df = current_row[, "df"], dispstr = "un"
  )

  # 指定边际分布的参数。对每个FAC_NAME 找到其在合并参数总表中的参数列。
  # 最后返回的list 每个对象都命名为"dp"
  margin_param_list <- lapply(FAC_NAMES, FUN = function(name) {
    current_row[, grep(name, colnames(current_row))]
  })
  names(margin_param_list) <- rep("dp", times = length(FAC_NAMES))

  # 通过t_copula 和skew t 边际分布生成多维分布，从而得到随机数。
  # 返回的结果起名为FAC_NAMES 
  row_mvdc <- mvdc(copula = t_cop_row, margins = c("st", "st", "st"), paramMargins = margin_param_list)
  set.seed(seed = seed)
  multidis_random_nums <- rMvdc(n = n, mvdc = row_mvdc)
  colnames(multidis_random_nums) <- FAC_NAMES

  return(multidis_random_nums)
}
head(randomNum_tCop_stMargin(merged_params[1, ], seed = 101))


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
