# ===================================================================
# 使用rmgarch 包进行garch-copula 模型进行拟合的测试脚本。
# 这里是对于整个时间框架内的数据进行拟合，分别实验了不同的copula 指定
# ===================================================================
library(rmgarch)
library(dplyr)
library(xts)
library(readr)

week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))
year_end_idx <- endpoints(week_fac, on = "year")

garch_type <- "eGARCH"
var_mdl <- list(model = garch_type, garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")

fit_garch_copula <- function(uarchSpec, fac_data, copula_type, dcc, asymm = FALSE) {
  # 输入
  # archSpec: 单变量arch 指定; fac_data: 因子数据框; 
  # copula_type: str c("mvnorm", "mvt"); dcc: Bool 是否为动态copula
  cop_spec <- 
    cgarchspec(uspec = multispec(replicate(3, uarchSpec)), 
               asymmetric = asymm,
               distribution.model = list(copula = copula_type, method = "ML", time.varying = dcc, transformation = "parametric"))
  cop_fit <- cgarchfit(cop_spec, data = fac_data)
  return(cop_fit)
  
  # 抽出与相关性有关的参数值
  # coefs_vec <- coef(cop_fit)
  # joint_para <- coefs_vec[grep("^\\[Joi.+", names(coefs_vec))]
  # 
  # # 输出极大似然函数值
  # likelihood_value <- likelihood(cop_fit)
  # 
  # return(c(joint_para, likelihood = likelihood_value))
}

static_tCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvt", dcc = FALSE)
static_normCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvnorm", dcc = FALSE)
dcc_tCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvt", dcc = TRUE)
dcc_normCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvnorm", dcc = TRUE)
adcc_tCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvt", dcc = TRUE, asymm = TRUE)
adcc_normCop <- fit_garch_copula(uarchSpec = archSpec, week_fac, copula_type = "mvnorm", dcc = TRUE, asymm = TRUE)

garch_cop_fits <-  list(static_t = static_tCop, static_norm = static_normCop, 
                    dcc_t = dcc_tCop, dcc_norm = dcc_normCop,
                    adcc_t = adcc_tCop, adcc_norm = adcc_normCop)
# cop_params <- lapply(cop_params, round, 4)



# =======================================================================================================
# 临时测试
# =======================================================================================================
cop_spec <- 
  cgarchspec(uspec = multispec(replicate(3, archSpec)), 
             distribution.model = list(copula = "mvt", method = "ML", time.varying = F, transformation = "parametric"))
cop_fit <- cgarchfit(cop_spec, data = week_fac)









