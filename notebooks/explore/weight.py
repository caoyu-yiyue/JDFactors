# %%
import pandas as pd
from scipy.optimize import minimize, Bounds
import numpy as np
from numba import jit

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
DATE_LIST: list = list(random_num['date'].drop_duplicates())
FAC_NUM = random_num.shape[1]

# %%
# 通过现在所在的date 找到下一个date，并在rf 数据框中找到rf 数据
current_date = first_group['date'][0]
nxt_date_idx = DATE_LIST.index(current_date) + 1
nxt_date = DATE_LIST[nxt_date_idx]
rf = rf_df.loc[nxt_date, 'Nrrwkdt']


# %%
@jit(nopython=True)
def expect_value(weight: np.array, rf: float, gamma: int,
                 fac_ndarr: np.ndarray):
    """
    计算 CRRA 期望收益

    Parameter:
    ----------
    weight:
        np.array
        权重序列
    rf:
        float
        下一期的无风险收益
    gamma:
        CRRA 中的gamma 值
    fac_df:
        np.ndarray
        只保存因子的ndArray，不包含时间列

    Return:
    -------
        float
        期望收益的相反数，因为要用到 minimize 函数中最大化期望
    """

    weight = np.asarray(weight)
    # 计算一个收益序列
    utility_array: np.array = (1 + rf + fac_ndarr.dot(weight))**(1 - gamma) / (
        1 - gamma)

    # 返回收益列序列的均值的相反数，从而最大化期望收益。
    expected_value = utility_array.mean()
    return -expected_value


# %%
def jacob(weight: np.array, rf, gamma, fac_array: np.ndarray):
    """
    针对目标函数计算梯度向量

    Parameter:
    ----------
    weight:
        np.array
        权重array
    rf:
        float
        下一期的无风险收益
    gamma:
        float
        CRRA 中的gamma
    fac_df:
            np.ndarray
            只保存因子的ndArray，不包含时间列

    Result:
    -------
    np.array
        表示梯度向量的np.array
    """

    weight = np.asarray(weight)
    # 计算每个导数的共同部分，命名为gradient_core
    gradient_core: np.array = (1 + rf + fac_array.dot(weight))**(-gamma)
    # 在每点（随机数）上，每个梯度分量，是共同部分乘以相应的因子值。
    gradient_array: np.ndarray = fac_array * gradient_core[:, np.newaxis]

    # 最终在每列上求平均，得到目标函数的梯度向量。
    return gradient_array.mean(axis=0)


# %%
def hess(weight: np.array, rf, gamma, fac_array: np.ndarray):
    """
    针对目标函数计算梯度向量

    Parameter:
    ----------
    weight:
        np.array
        权重array
    rf:
        float
        下一期的无风险收益
    gamma:
        float
        CRRA 中的gamma
    fac_df:
            np.ndarray
            只保存因子的ndArray，不包含时间列

    Result:
    -------
    np.ndarray
        一个2d array，是对 weight 向量的海塞矩阵，
    """

    weight = np.asarray(weight)
    # 首先计算每个海塞矩阵分量需要的共同部分，每一个随机数对应一个，形成一个array
    hass_core = -gamma * (1 + rf + fac_array.dot(weight))**(-gamma - 1)

    # 计算海塞矩阵分量的后半部分，是fac_array 每一行的向量外积。
    # 为了使每一行的计算自动归为一个ndarray，将fac_array 升到 3 维，
    # 乘号前面的分量是一个列向量，后面的分量为行向量，最终在每个分量处对应相乘，得到向量的外积。
    outer_product_array = fac_array[:, :, np.newaxis] * fac_array[:, np.
                                                                  newaxis, :]
    # 每个hass_core 分量乘对应的外积矩阵分量。这里将hass_core 升到 3 维相乘即可。
    hass_matr_array: np.ndarray = outer_product_array * hass_core[:, np.
                                                                  newaxis, np.
                                                                  newaxis]

    # 最终计算每个随机数形成的海塞矩阵的平均值，得到原 objective 函数的海塞矩阵
    return hass_matr_array.mean(axis=0)


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
bnd = Bounds([0] * (FAC_NUM - 1), [1] * (FAC_NUM - 1))

# %%
first_array = first_group.drop('date', axis=1).to_numpy()

# %%
sol = minimize(fun=expect_value,
               x0=[0., 0., 0.],
               args=(rf, 7, first_array),
               method='trust-constr',
               jac=jacob,
               hess=hess,
               constraints=con,
               bounds=bnd)
sol
