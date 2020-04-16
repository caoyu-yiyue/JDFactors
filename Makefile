.PHONY: all prepare_data garch_model

all: data/interim/fac_xts_Week.Rds data/interim/fac_xts_Week.Rds garch_model

clean:
	trash data/interim/*.Rds
	trash data/processed/*.Rds

# ================================= prepare data =================================== #
# 该脚本同时产生三个文件，这里选择其中一个作为target
data/interim/fac_xts_Week.Rds:
	Rscript --vanilla src/data/prepare_data.R

# ================================= garch model =================================== #
garch_model: data/interim/best_arma_ssdt_Week.Rds data/interim/multi_garch_mdl.Rds

data/interim/best_arma_ssdt_Week.Rds: data/interim/fac_xts_Week.Rds
	Rscript --vanilla src/process/best_arma.R --data_freq "Week" -o $@

data/interim/multi_garch_mdl.Rds: data/interim/best_arma_ssdt_Week.Rds
	Rscript --vanilla src/process/multi_garch_mdl.R $@
