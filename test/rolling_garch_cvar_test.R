suppressPackageStartupMessages({
  library(xts)
})

source("src/data/read_data.R")
source("src/process/multi_garch_mdl.R")

rolling_garch_list <- read_rolling_multigarchfit()
rolling_garch_list <- readRDS("data/interim/rolling_multigarch.Rds")


test_conv_cvar <- function(mgfit) {
  conv_flags <- .conver_for_multigarchfit(mgfit)
  cvar_flags <- .cvar_for_multigarchfit(mgfit)
  if (all(conv_flags) & all(cvar_flags)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


result <- sapply(rolling_garch_list, test_conv_cvar)
which(!result)
