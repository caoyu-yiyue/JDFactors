.PHONY: all best_weight

all: best_weight

seed:= 10001
nbins:= 100
gamma:= 7
mr:= 0.2
# half:= first
method:= DE
max_r:= 0.5
sum_1:= T
sum_1_flag= $(if $(filter T, $(sum_1)), --sum_1, --no_sum_1)

data_dir:= data/processed/ga$(gamma)_mr$(mr)_max$(max_r)_sum1$(sum_1)/
best_weight_fpath_first:= $(data_dir)best_weight_s$(seed)_nb$(nbins)_m$(method)_half1st.pickle
best_weight_fpath_second:= $(data_dir)best_weight_s$(seed)_nb$(nbins)_m$(method)_half2nd.pickle

$(best_weight_fpath_first): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --mr $(mr) --half first --method $(method) --max_r $(max_r) $(sum_1_flag) $@
$(best_weight_fpath_second): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --mr $(mr) --half second --method $(method) --max_r $(max_r) $(sum_1_flag) $@
best_weight_first: $(best_weight_fpath_first)
best_weight_second: $(best_weight_fpath_second)

# 执行如下命令会同时开两个python 进程，但是只会显示一个进度条
best_weight: 
	$(MAKE) best_weight_first &
	$(MAKE) best_weight_second



############################## 合并数据框 ###############################
concat_halfs:
	python3 /Users/caoyue/Codes/Python/JDFactors/src/process/combine_halfs.py
combine_by_mins: concat_halfs
	python3 /Users/caoyue/Codes/Python/JDFactors/src/process/combine_mins.py
	