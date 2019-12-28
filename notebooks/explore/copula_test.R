
library(rugarch)
library(copula)
library(readr)
library(dplyr)
library(tidyr)
library(xts)
library(sn)
library(foreach)
library(timetk)

# 读取数据并转换成xts
week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))
year_end_idx <- endpoints(week_fac, on = "year")


# garch 模型的指定
garch_type <- "eGARCH"
var_mdl <- list(model = garch_type, garchOrder = c(1, 1))
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
arch_roll_dist_path <- paste("data/interim/", garch_type, "_arch_roll_dist.Rda", sep = "")
save(arch_roll, arch_roll_dist_df, file = arch_roll_dist_path)
load(file = arch_roll_dist_path)

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
  
  ############## 对该段窗口进行 garchfit 并得到 standarize residuals ##################
  garch_fit_list <- lapply(d_window, FUN = function(col) {
    fit <- ugarchfit(archSpec, data = col, solver = "hybrid")
    return(fit)
  })
  # 获取到每个因子的standardize residuals
  resids_list <- lapply(garch_fit_list, FUN = function(fit) {
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

  ###################### 估计 copula 系数 ###############################3
  pob_var <- pobs(as.matrix(resids_xts))
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
merged_params_path <- paste("data/interim/", garch_type, "_merged_params.Rda", sep = "")
save(merged_params, file = merged_params_path)
load(merged_params_path)




########################################################################################
#########                          生成随机数                            ###############
########################################################################################
FAC_NAMES <- names(week_fac)
gen_mvdcs_tCop_stMargin <- function(current_row) {
  # 根据保存边际分布参数与copula 参数的行，生成联合分布的随机数
  # current_row: 保存参数的行。这里apply 会传入一个numeric vector
  # seed: 生成随机数的seed，默认为101。
  # n: 需要的随机数的个数，默认为100000。
  
  # 每一行单独进行操作，所以不需要index。否则后续会把index 认为是多的一列。
  current_row <- coredata(current_row)
  
  # 指定t_copula 对象
  t_cop_row <- tCopula(
    param = current_row[c("rho.1", "rho.2", "rho.3")],
    dim = length(FAC_NAMES), df = current_row["df"], dispstr = "un"
  )
  
  # 指定边际分布的参数。对每个FAC_NAME 找到其在合并参数总表中的参数列。
  # 最后返回的结构为list-list$dp，每个名为dp 的个体下是一个vector，保存了四个skew-t 参数。
  margin_param_list <- lapply(FAC_NAMES, FUN = function(name) {
    row_single_fac <- current_row[grep(name, names(current_row))]
    list_single_fac <- list(dp = row_single_fac)
    return(list_single_fac)
  })
  
  # 通过t_copula 和skew t 边际分布生成多维分布，从而得到随机数。
  # 返回的结果起名为FAC_NAMES
  row_mvdc <- mvdc(copula = t_cop_row, margins = c("st", "st", "st"), paramMargins = margin_param_list)
  
  return(row_mvdc)
}

mvdc_list <- apply(merged_params, gen_mvdcs_tCop_stMargin, MARGIN = 1)

# test_random_num <- randomNum_tCop_stMargin(merged_params[1, ], seed = 101)
library(doSNOW)
gen_multi_dist_random <- function(mvdcs, seed, n, save_file = FALSE) {
  row_n <- length(mvdcs)
  pb <- txtProgressBar(max = row_n, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # 使用并行循环计算随机数（按行数循环）
  
  multi_dist_random_num <- foreach(i = 1:row_n, .combine = "rbind", .packages = c("xts", "copula", "sn"), 
                                   .export = c("FAC_NAMES"),
                                   .options.snow = opts) %dopar% {
    set.seed(n)
    multidis_random_nums <- rMvdc(n = n, mvdc = mvdcs[[i]])
    colnames(multidis_random_nums) <- FAC_NAMES
    
    return(multidis_random_nums)
  }
  close(pb)

  # 加入date 列，成为一个data.frame 然后输出
  random_num_df <- as.data.frame(multi_dist_random_num)
  date <- rep(names(mvdcs), each = n)
  result_df <- cbind(date, random_num_df)
  
  if(save_file) {
    random_num_path <- paste("data/interim/", garch_type, "_random_num_", seed, ".Rds", sep = "")
    saveRDS(result_df, file = random_num_path)
  }
  return(result_df)
}


# 随机数的生成，浪费时间！！！！！！
seeds = c(101, 304, 10001, 100001, 1000001)
# random_num_result_1 <- gen_multi_dist_random(param_df = merged_params, seed = seeds[1],
#                                            n = 10000, save_file = TRUE)

# 开一个结果df，同时给列命名
result_rands <- data.frame(matrix(nrow = 0, ncol = length(FAC_NAMES)))
colnames(result_rands) <- FAC_NAMES

inte_num <- 9 # 总共需要计算多少
now = 0 # 目前已经计算了多少
# 从某个seed 出发，加上一个总数，在inte_num 的范围内生出新seed 做循环
for(s in (seeds[3] + now):(seeds[3] + inte_num)) {
  # 打印当前的循环个数，以及完成的百分比
  counter <- s - seeds[3]
  print(sprintf("NO. %d, %.2f%%", counter + 1, (counter * 100) / inte_num))
 
  # 开一个cluster，然后并行在mvdc_list 上循环生成随机数1000 个。
  cls <- makeCluster(4, type = "SOCK")
  registerDoSNOW(cls)
  part_rands <- gen_multi_dist_random(mvdcs = mvdc_list, seed = s,
                                       n = 1000, save_file = FALSE)
  stopCluster(cls)
  gc()
  # 与结果df 合并，从而将结果输出出来。
  result_rands <- rbind(result_rands, part_rands)
}
saveRDS(result_rands, file = paste("data/interim/", garch_type, "_random_num_", seeds[3], ".Rds", sep = ""))


# 合并分开的随机数结果，并保存为一个feather 结果 #
path_patter <- paste(garch_type, "_random_num_\\d+.Rds", sep = "")
random_num_paths <- list.files("data/interim", pattern = path_patter)
random_num_list <- lapply(random_num_paths, FUN = function(path) {
  full_path <- paste("data/interim/", path, sep = "")
  result <- readRDS(full_path)
  return(result)
})

# date 列保存为Date 类
random_num_list <- lapply(random_num_list, FUN = function(df) { 
  df$date <- lubridate::ymd(df$date)
  return(df)
})
binded_random_df <- Reduce(f = function(top, bott) {rbind(top, bott)}, x = random_num_list)

# 保存为.feather
library(feather)
binded_path <- paste("data/interim/", garch_type, "_random_num_all.feather", sep = "")
write_feather(x = binded_random_df, path = binded_path)













# # 使用生成的随机数最大化权重值
# library(PortfolioAnalytics)
# CRRA <- function(R, weights, lambda, sigma, m3, m4) {
#   weights <- matrix(weights, ncol = 1)
#   M2.w <- t(weights) %*% sigma %*% weights
#   M3.w <- t(weights) %*% m3 %*% (weights %x% weights)
#   M4.w <- t(weights) %*% m4 %*% (weights %x% weights %x% weights)
#   term1 <- (1 / 2) * lambda * M2.w
#   term2 <- (1 / 6) * lambda * (lambda + 1) * M3.w
#   term3 <- (1 / 24) * lambda * (lambda + 1) * (lambda + 2) * M4.w
#   out <- -term1 + term2 - term3
#   out
# }
# 
# crra.moments <- function(R, ...) {
#   out <- list() + out$sigma <- cov(R)
#   out$m3 <- PerformanceAnalytics:::M3.MM(R)
#   out$m4 <- PerformanceAnalytics:::M4.MM(R)
#   out
# }
# 
# crra_port <- portfolio.spec(assets = colnames(test_random_num))
# crra_port <- add.constraint(crra_port, type = "full_investment")
# crra_port <- add.objective(portfolio = crra_port, type = "return", name = "CRRA",
#                            arguments = list(lambda = 10))
# opt_crra <- optimize.portfolio(R = test_random_num, portfolio = crra_port,
#                                momentFUN = "crra.moments", 
#                                optimize_method = "DEoptim")

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
