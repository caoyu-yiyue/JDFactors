source("src/result/result_tables.R")
source("src/data/read_data.R")

port_list <- read_port_ret()
weights_list <- read_opt_weights(which = "all")
facs_xts <- read_fac_xts()
rf_xts <- read_rf_xts()

result_table_sum1 <- result_table_main(
  port_ret_list = port_list$sum1,
  weights_list = weights_list$sum1, facs_xts = facs_xts
)

result_table_nosum1 <- result_table_main(
  port_ret_list = port_list$no_sum1,
  weights_list = weights_list$no_sum1, facs_xts = facs_xts, rf_xts = rf_xts
)
