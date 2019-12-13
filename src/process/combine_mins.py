# %%
import pandas as pd
from glob import glob
from functools import reduce

# %%
weights_root_dir = 'data/processed/'
# 根目录下的子目录列表
weight_sub_dirs = glob(weights_root_dir + 'ga*/')


def compare_combine(df_left, df_right):
    # 生成一列bool 值，指示左侧是否小于右侧的数据框
    left_min_T = df_left['func_v'] < df_right['func_v']

    # 左侧小的，留下左侧数据框中的所有行；右侧小的，留下右侧小的行
    left_remain = df_left[left_min_T]
    right_remain = df_right[~left_min_T]

    # 左右留下的合并到一起
    remain_df: pd.DataFrame = pd.concat([left_remain, right_remain])
    remain_df.sort_index(inplace=True)
    return remain_df


# 对所有的子目录进行循环
for sub_dir in weight_sub_dirs:
    # 找到其中的all 结尾的文件
    files_in_sub = glob(sub_dir + '*all.pickle')
    dfs_in_sub = map(pd.read_pickle, files_in_sub)
    # print(list(dfs_in_sub)[0])
    combined_df: pd.DataFrame = reduce(compare_combine, dfs_in_sub)
    combined_df.to_pickle(sub_dir + 'combined_result.pickle')
