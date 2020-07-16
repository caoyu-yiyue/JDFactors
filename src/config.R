#############################################################################
# 一些项目中需要使用的常量
#############################################################################

ROLLING_STEP <- c(Week = 12, Day = 60, Month = 12)
# 对于Month，这里使用的是行号，而非年数，用于精细控制起始行
IN_SAMPLE_YEARS <- c(Week = 8, Day = 8, Month = 100)

ROLLING_ARMA_ORDERS <- matrix(rep(3, 10), nrow = 2)
ROLLING_GARCH_ORDERS <- matrix(rep(1, 10), nrow = 2)
