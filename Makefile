.PHONY: all prepare_data

all: prepare_data data/interim/best_arma_ssdt_Week.Rds

prepare_data:
	Rscript --vanilla src/data/prepare_data.R

data/interim/best_arma_ssdt_Week.Rds: prepare_data
	Rscript --vanilla src/process/best_arma.R --data_freq "Week" -o $@
