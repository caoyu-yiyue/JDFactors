.PHONY: all best_weight

all: best_weight

seed:= 10001
nbins:= 100
gamma:= 7
# half:= first
method:= DE
max_r:= 0.5

best_weight_fpath_first:= data/processed/best_weight_s$(seed)_nb$(nbins)_ga$(gamma)_half1st_m$(method)_max$(max_r).pickle
best_weight_fpath_second:= data/processed/best_weight_s$(seed)_nb$(nbins)_ga$(gamma)_half2nd_m$(method)_max$(max_r).pickle
$(best_weight_fpath_first): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --half first --method $(method) --max_r $(max_r) $@
$(best_weight_fpath_second): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --half second --method $(method) --max_r $(max_r) $@
best_weight_first: $(best_weight_fpath_first)
best_weight_second: $(best_weight_fpath_second)

# 执行如下命令会同时开两个python 进程，但是只会显示一个进度条
best_weight: 
	$(MAKE) best_weight_first &
	$(MAKE) best_weight_second
