# %%
import pandas as pd
from glob import glob
# import os

# %%
weights_root_dir = 'data/processed/'
# processed 目录下的子目录列表
weights_sub_dirs: list = glob(weights_root_dir + '*/')


def get_unique_prefix(file_name: str):
    if 'half' in file_name:
        return file_name.split('_half')[0]
    else:
        return None


# %%
# 对子目录列别下的每个目录循环
for sub_dir in weights_sub_dirs:
    # 每个子目录中的所有文件列表
    files_in_sub = glob(sub_dir + '*')

    # 所有文件砍掉 half1st.pickle 和half2nd.pickle，所有的前缀去重形成的列表
    file_prefixs: list = list(set(map(get_unique_prefix, files_in_sub)))
    # 如果前缀列表中有None 则去掉
    if None in file_prefixs:
        file_prefixs.remove(None)

    # 对前缀列表中的每个前缀循环
    for file_pre in file_prefixs:
        # 每个前缀加上尾巴，读取文件然后合并、保存整个数据框的文件
        file_path_1st = file_pre + '_half1st.pickle'
        file_path_2nd = file_pre + '_half2nd.pickle'

        df_1st = pd.read_pickle(file_path_1st)
        df_2nd = pd.read_pickle(file_path_2nd)
        df_all: pd.DataFrame = pd.concat([df_1st, df_2nd])

        df_all.to_pickle(file_pre + '_all.pickle')
        # os.remove(file_path_1st)
        # os.remove(file_path_2nd)
