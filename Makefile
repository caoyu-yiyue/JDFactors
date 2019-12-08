.PHONY: all best_weight

all: best_weight

seed:= 10001
nbins:= 2000
gamma:= 7

best_weight_fpath:= data/interim/best_weight_s$(seed)_nb$(nbins)_ga$(gamma).pickle
$(best_weight_fpath): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma)
best_weight: $(best_weight_fpath)
