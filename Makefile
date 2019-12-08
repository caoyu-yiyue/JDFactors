.PHONY: all best_weight

all: best_weight

seed:= 10001
nbins:= 100
gamma:= 7
half:= first

best_weight_fpath:= data/interim/best_weight_s$(seed)_nb$(nbins)_ga$(gamma)_half$(half).pickle
$(best_weight_fpath): data/interim/eGARCH_random_num_all.feather | data/interim/eGARCH_random_num_all.feather
	python3 -O src/optimize_weight.py --seed $(seed) --nbins $(nbins) --gamma $(gamma) --half $(half)
best_weight: $(best_weight_fpath)
