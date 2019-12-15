library(rugarch)
library(readr)
library(dplyr)
library(xts)
library(feather)
library(sn)
library(copula)

week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))

# 建立garch 模型
garch_type <- "eGARCH"
var_mdl <- list(model = garch_type, garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")

# 获取到每个因子的standardize residuals
resids_list <- lapply(week_fac, FUN = function(col) {
  fit <- ugarchfit(archSpec, data = col)
  return(residuals(fit, standardize = TRUE))
})

# 每个fac 加名字
fac_names <- names(resids_list)
for (fac in fac_names) {
  colnames(resids_list[[fac]]) <- fac
}
# 合并每个因子的residual 成为xts 对象
resids_xts <- Reduce(f = function(left_xts, right_xts) {
  merge.xts(left_xts, right_xts)
}, x = resids_list)

######################## copula 模型拟合出随机数 ##########################
# 估计分布
start_pnt <- nrow(week_fac) - 2

# 使用所有的数据，估计出最后一天的分布参数
garch_density <- lapply(week_fac, function(col) {
  fit <- ugarchroll(archSpec,
    data = col, n.start = start_pnt,
    refit.window = "recursive", refit.every = 1,
    calculate.VaR = FALSE, keep.coef = FALSE
  )
  density_df <- as.data.frame(fit, which = "density")
  density_param <- density_df[nrow(density_df), 1:4]
  return(density_param)
})

# 将cp 参数转换为dp 参数
cpDf2dpDf <- function(df) {
  # 对df 中的每一行进行cp2dp，最后返回一个xts 对象
  dp_matrix <- t(apply(df, 1, FUN = cp2dp, family = "ST"))

  # 确保index 为Date 对象，不包括时间
  idx <- as.Date(rownames(dp_matrix))
  return(xts(dp_matrix, order.by = idx))
}
dp_density_params <- lapply(garch_density, cpDf2dpDf)
margin_params <- lapply(dp_density_params, FUN = function(param_xts) {
  coredata(param_xts)
})

# 估计copula 参数
pob_vars <- pobs(as.matrix(week_fac))
d <- ncol(week_fac)
cop_mdl <- tCopula(dim = d, dispstr = "un")
cop_fit <- fitCopula(cop_mdl, pob_vars, method = "ml")

rho <- coef(cop_fit)[1:3]
df <- coef(cop_fit)[4]

# 估计出联合分布
tcop_stMargin_dist <- mvdc(
  copula = tCopula(rho, dim = d, df = df, dispstr = "un"), margins = c("st", "st", "st"),
  paramMargins = list(
    list(dp = margin_params[["premium"]]), list(dp = margin_params[["smb"]]),
    list(dp = margin_params[["hml"]])
  )
)
# 生成和week_fac 行数一样多的随机数
set.seed(10001)
tcop_stMargin_rand_nums <- rMvdc(nrow(week_fac), tcop_stMargin_dist)
colnames(tcop_stMargin_rand_nums) <- colnames(week_fac)


################################ norm copula ################################
norm_cop <- normalCopula(dim = d, dispstr = "un")
norm_cop_fit <- fitCopula(norm_cop, pob_vars, method = "ml")

norm_rhos <- coef(norm_cop_fit)
norm_cop_stMargin_dist <- mvdc(
  copula = tCopula(norm_rhos, dim = d, dispstr = "un"), margins = c("st", "st", "st"),
  paramMargins = list(
    list(dp = margin_params[["premium"]]), list(dp = margin_params[["smb"]]),
    list(dp = margin_params[["hml"]])
  )
)

set.seed(10001)
normcop_stMargin_rand_nums <- rMvdc(10000, norm_cop_stMargin_dist)
colnames(normcop_stMargin_rand_nums) <- colnames(week_fac)

############################## 保存结果 #####################################
# resids_xts 转化为df 并保存为feather
resids_df <- data.frame(date = index(resids_xts), coredata(resids_xts))
write_feather(resids_df, path = "data/interim/garch_residuals.feather")

# 保存通过garch 和copula 模型生成的随机数
write_feather(as.data.frame(tcop_stMargin_rand_nums), path = "data/interim/garch_tcop_stMargin_randoms.feather")
write_feather(as.data.frame(normcop_stMargin_rand_nums), path = "data/interim/garch_normcop_stMargin_randoms.feather")
