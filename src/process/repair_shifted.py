# %%
import pandas as pd
from glob import glob
import re
import os
from datetime import timedelta
from numpy import nan
from src.process import optimize_weight as ow
from mystic.symbolic import (generate_conditions, generate_constraint,
                             generate_penalty, generate_solvers, simplify)

# %%
rand_num_df: pd.DataFrame = ow.read_simul_random_df()
rf_df: pd.DataFrame = ow.read_rf_df()
merge_df = pd.DataFrame = ow.join_rf_df(rand_num_df, rf_df, rf_type='week')
FAC_NAME = rand_num_df.columns.to_list()
FAC_NUM = rand_num_df.shape[1]

# %%
shifted_files = glob('data/processed/shifted/*/*')

for shifted_f in shifted_files:
    # %%
    m = re.search(
        (r'ga(?P<gamma>\d)_mr(?P<mr>.+)_max(?P<max_r>\d)_sum1(?P<sum1>[T|F])'
         r'/.*_s(?P<seed>\d+)_nb(?P<nbins>\d+)_m(?P<method>\w+)_half'),
        shifted_f)
    gamma = int(m.group('gamma'))
    mr = float(m.group('mr'))
    max_r = float(m.group('max_r'))
    sum1 = True if m.group('sum1') == 'T' else False
    seed = int(m.group('seed'))
    nbins = int(m.group('nbins'))
    method = m.group('method')

    # %%
    # 末尾加一行nan，然后整体向下移动一行
    ini_df = pd.read_pickle(shifted_f)
    last_idx = ini_df.index[-1] + timedelta(days=7)
    cols = ini_df.columns
    last_row = pd.Series(data=[nan] * len(cols), index=cols, name=last_idx)

    add_last_df = ini_df.append(last_row)
    first_row_blank: pd.DataFrame = add_last_df.shift(1)

    # %%
    first_idx = ini_df.index[0]
    first_group_ran_nums: pd.DataFrame = merge_df.loc[first_idx, :]

    equation = ow.equation_str(mr=mr, sum_1=sum1)
    pf = generate_penalty(generate_conditions(equation), k=1e20)
    cf = generate_constraint(generate_solvers(simplify(equation)))

    first_row_weight: pd.Series = ow.opti_fun(df=first_group_ran_nums,
                                              gamma=gamma,
                                              max_r=max_r,
                                              seed=seed,
                                              nbins=nbins,
                                              method=method,
                                              constraint=cf,
                                              penalty=pf,
                                              FAC_NAME=FAC_NAME,
                                              FAC_NUM=FAC_NUM)

    first_row_blank.loc[first_idx] = first_row_weight
    repaired_df = first_row_blank

    save_path = shifted_f.replace('shifted/', '')
    save_dir = os.path.split(save_path)[0]
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    repaired_df.to_pickle(save_path)
