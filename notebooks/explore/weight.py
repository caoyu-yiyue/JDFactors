# %%
import pandas as pd
from scipy.optimize import minimize

# %%
random_num: pd.DataFrame = pd.read_feather('data/interim/random_num.feather')
random_num['date'] = pd.to_datetime(random_num['date'], format='%Y-%m-%d')
grouped = random_num.groupby("date")
first_group: pd.DataFrame = grouped.get_group(('2003-01-03'))

# %%
rf_df: pd.DataFrame = pd.read_csv('data/raw/csvFiles/TRD_Nrrate.csv',
                                  sep='\t',
                                  parse_dates=['Clsdt'],
                                  index_col='Clsdt',
                                  encoding='utf-16')

# %%
date_list: list = list(random_num['date'].drop_duplicates())


# %%
def expect_value(weight, gamma, fac_df):
    """
    计算 CRRA 期望收益

    Parameter:
    ----------
    weight:
        list
        权重序列
    gamma:
        CRRA 中的gamma 值
    fac_df:
        带有factor 收益随机数的DataFrame，其中带一列date 表示日期。

    Return:
    -------
        float
        期望收益的相反数，因为要用到 minimize 函数中最大化期望
    """
    # 通过现在所在的date 找到下一个date，并在rf 数据框中找到rf 数据
    current_date = first_group['date'][0]
    nxt_date_idx = date_list.index(current_date) + 1
    nxt_date = date_list[nxt_date_idx]
    rf = rf_df.loc[nxt_date, 'Nrrwkdt']

    # 去掉date 列，只剩下因子列
    facs_data: pd.DataFrame = first_group.drop('date', axis=1)

    # 计算一个收益序列
    return_series: pd.Series = ((facs_data.mul(weight, axis=1).sum(axis=1) +
                                 rf + 1)**(1 - gamma)) / (1 - gamma)

    # 返回收益列序列的均值的相反数，从而最大化期望收益。
    expected_value = return_series.mean()
    return -expected_value


# %%
# 限制条件
def weight_constrain(weight, mr):
    """
    对权重的和施加限制条件，
    Parameters:
    -----------
    weight:
        list
        权重序列
    mr:
        float
        杠杆率
    """
    return (1 / mr) - weight[0] - (2 * sum(weight[1:]))


con = {'type': 'ineq', 'fun': weight_constrain, 'args': (0.2, )}

# %%
bnd = ([0, 1], ) * (first_group.shape[1] - 1)

# %%
sol = minimize(fun=expect_value,
               method='SLSQP',
               x0=[0., 0., 0.],
               constraints=con,
               args=(7, first_group),
               bounds=bnd)
sol
