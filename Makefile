.PHONY: all best_weight

all: best_weight

seed:= 10001
nbins:= 100
gamma:= 7
# half:= first
method:= DE

best_weight_fpath_first:= data/processed/best_weight_s$(seed)_nb$(nbins)_ga$(gamma)_half1st_m$(method).pickle
best_weight_fpath_second:= data/processed/best_weight_s$(seed)_nb$(nbins)_ga$(gamma)_half2nd_m$(method).pickle
$(best_weight_fpath_first): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --half first --method $(method) $@
$(best_weight_fpath_second): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --half second --method $(method) $@
best_weight_first: $(best_weight_fpath_first)
best_weight_second: $(best_weight_fpath_second)
best_weight: best_weight_first best_weight_second
