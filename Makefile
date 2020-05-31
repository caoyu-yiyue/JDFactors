.PHONY: all prepare_data garch_model

all: data/interim/facs_xts.Rds data/interim/rf_xts.Rds garch_model rolling_part result

clean:
	trash data/interim/*.Rds
	trash data/processed/*.Rds

# ================================= prepare data =================================== #
# 该脚本同时产生三个文件，这里选择其中一个作为target
data/interim/facs_xts.Rds:
	Rscript --vanilla src/data/prepare_data.R

data/interim/rf_xts.Rds:
	Rscript --vanilla src/data/rf_data.R

# ================================= garch model =================================== #
garch_model: data/interim/best_arma_ssdt_Week.Rds data/interim/multi_garch_mdl.Rds \
	data/processed/all_cops.Rds data/interim/opt_weights.Rds

data/interim/best_arma_ssdt_Week.Rds: data/interim/facs_xts.Rds
	Rscript --vanilla src/process/best_arma.R --data_freq "Week" -o $@

data/interim/multi_garch_mdl.Rds: data/interim/best_arma_ssdt_Week.Rds
	Rscript --vanilla src/process/multi_garch_mdl.R $@

data/processed/all_cops.Rds: data/interim/multi_garch_mdl.Rds
	Rscript --vanilla src/process/copula_mdl.R $@

# ================================= rolling fit =================================== #
data_freq:= Week

rolling_part: data/interim/rolling_multigarch_$(data_freq).Rds \
	data/interim/rolling_cop_rcov_$(data_freq).Rds data/interim/fixed_cor2cov_$(data_freq).Rds \
	data/interim/rolling_mean_$(data_freq).Rds data/interim/opt_weights_$(data_freq).Rds \
	data/processed/port_ret_$(data_freq).Rds

data/interim/rolling_multigarch_$(data_freq).Rds: data/interim/facs_xts.Rds
	Rscript --vanilla src/process/rolling_multigarch.R -f $(data_freq) -o $@

data/interim/rolling_cop_rcov_$(data_freq).Rds: data/interim/rolling_multigarch_$(data_freq).Rds
	Rscript --vanilla src/process/rolling_copula.R -f $(data_freq) -o $@

data/interim/fixed_cor2cov_$(data_freq).Rds: data/interim/rolling_multigarch_$(data_freq).Rds
	Rscript --vanilla src/process/fixed_cor2cov.R -f $(data_freq) -o $@

data/interim/rolling_mean_$(data_freq).Rds: data/interim/facs_xts.Rds
	Rscript --vanilla src/process/rolling_mean.R -f $(data_freq) -o $@

data/interim/opt_weights_$(data_freq).Rds: data/interim/rolling_cop_rcov_$(data_freq).Rds \
	data/interim/fixed_cor2cov_$(data_freq).Rds data/interim/rolling_mean_$(data_freq).Rds
	Rscript --vanilla src/process/weight_optimize.R -f $(data_freq) -o $@

# ================================ inverst result =================================== #
data/processed/port_ret_$(data_freq).Rds: data/interim/opt_weights_$(data_freq).Rds
	Rscript --vanilla src/result/port_ret.R -f $(data_freq) -o $@
