library(rugarch)
library(rmgarch)
library(xts)
library(readr)
library(dplyr)

week_raw <- read_csv("data/raw/csvFiles/week_readed.csv")
week_use <- filter(week_raw, trdWeek > "2000-01-01")
week_fac <- as.xts(week_use[, 2:4], order.by = as.Date(week_use$trdWeek))
year_end_idx <- endpoints(week_fac, on = "year")

garch_type <- "eGARCH"
var_mdl <- list(model = garch_type, garchOrder = c(1, 1))
mean_mdl <- list(armaOrder = c(3, 0))
archSpec <- ugarchspec(variance.model = var_mdl, mean.model = mean_mdl, distribution.model = "sstd")

multi_spec <- multispec(replicate(3, archSpec))
multi_fit <- multifit(multispec = multi_spec, data = week_fac)

cop_spec <- cgarchspec(uspec = multi_spec, 
               distribution.model = list(copula = "mvt", method = "ML", time.varying = FALSE, transformation = "parametric"))
cop_fit <- cgarchfit(cop_spec, data = week_fac, fit = multi_fit)


  
  
  
  
  
  