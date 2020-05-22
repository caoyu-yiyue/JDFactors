suppressPackageStartupMessages({
  library(xts)
})

source("src/process/multi_garch_mdl.R")


test_conv_cvar <- function(mgfit) {
  conv_flags <- .conver_for_multigarchfit(mgfit)
  cvar_flags <- .cvar_for_multigarchfit(mgfit)
  if (all(conv_flags) & all(cvar_flags)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
