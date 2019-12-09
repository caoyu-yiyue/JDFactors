# %%
import pandas as pd

# %%
first_half = pd.read_pickle(
    'data/interim/best_weight_s10001_nb100_ga7_halffirst.pickle')
second_half = pd.read_pickle(
    'data/interim/best_weight_s10001_nb100_ga7_halfsecond.pickle')
result_df = pd.concat([first_half, second_half])

# %%
facs_ret = pd.read_csv('data/raw/csvFiles/week_readed.csv',
                       index_col='trdWeek',
                       parse_dates=['trdWeek'])
facs_ret.rename_axis(index={'trdWeek': 'date'}, inplace=True)
facs_ret = facs_ret.reindex(result_df.index)
FAC_NAMES = facs_ret.columns.to_list()

# %%
weights_df: pd.DataFrame = result_df[FAC_NAMES]
ret_series: pd.Series = weights_df.mul(facs_ret,
                                       axis='index').sum(axis='columns')

(ret_series + 1).prod() - 1
