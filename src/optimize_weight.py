# %%
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import numpy as np
import pandas as pd
import click
from mystic.solvers import LatticeSolver, PowellDirectionalSolver
from mystic.symbolic import (generate_conditions, generate_constraint,
                             generate_penalty, generate_solvers, simplify)
from mystic.termination import ChangeOverGeneration as COG
from mystic.tools import random_seed
from mystic.pools import SerialPool as Pool
# from distributed import Client


def read_simul_random_df(file='data/interim/eGARCH_random_num_all.feather'):
    """
    读取从r 语言模拟生成的多维随机数
    Parameters:
    ------------
    file: str
        读取文件的路径。这里是一个.feather 文件

    Return:
    -------
    pd.DataFrame
        读取到的数据框
    """
    random_num: pd.DataFrame = pd.read_feather(file)
    random_num['date'] = pd.to_datetime(random_num['date'], format='%Y-%m-%d')
    random_num.set_index('date', inplace=True)
    random_num.sort_index(inplace=True)
    return random_num


def read_rf_df(file='data/raw/csvFiles/TRD_Nrrate.csv'):
    """
    读取无风险收益数据，以时间（日）为index，三列分别为日度、周度、月度化的无风险收益
    Parameters:
    ----------
    file: str
        一个保存无风险收益率的.csv 文件

    Return:
    pd.DataFrame
        读取到的无风险收益数据框
    """
    rf_df: pd.DataFrame = pd.read_csv(file,
                                      sep='\t',
                                      header=0,
                                      encoding='utf-16le',
                                      usecols=[0, 2, 3, 4],
                                      names=['date', 'day', 'week', 'month'],
                                      parse_dates=['date'],
                                      index_col='date')
    rf_df.sort_index(inplace=True)
    return rf_df


def join_rf_df(random_num_df: pd.DataFrame, rf_df: pd.DataFrame, rf_type: str):
    """
    合并随机数数据框和一列无风险收益数据
    Parameters:
    -----------
    random_num_df:
        pd.DataFrame
        多维随机数数据框
    rf_df:
        pd.DataFrame
        无风险收益数据框
    rf_type:
        str, one list of ('day', 'week', 'month')

    Return:
    -------
    pd.DataFrame
        合并了一列无风险收益的数据框
    """
    # 将第二天的rf 提前一天，因为前一天需要使用后一天的rf 计算。
    shifted_rf: pd.Series = rf_df[rf_type].shift(-1)
    joined_df: pd.DataFrame = random_num_df.join(shifted_rf, how='inner')
    renamed_rf_col_df: pd.DataFrame = joined_df.rename(columns={rf_type: 'rf'})
    return renamed_rf_col_df


def _expect_value(weight, rf, gamma, fac_ndarr):
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


def equation_str(mr):
    equations = """
    x0 + x1 + x2 == 1
    x0 + 2*(x1 + x2) - (1. / {}) <= 0
    """.format(mr)

    return equations


def opti_latti_pow(df: pd.DataFrame, nbins, gamma, constraint, penalty, seed):
    # print(df.name)
    core_fac = df.drop(columns='rf')
    fac_array = core_fac.to_numpy()
    # fac_names = core_fac.columns
    rf = df['rf'][0]
    random_seed(s=seed)

    def helper_cost(weight):
        return _expect_value(weight, rf=rf, gamma=7, fac_ndarr=fac_array)

    solver = LatticeSolver(dim=FAC_NUM, nbins=nbins)
    solver.SetNestedSolver(PowellDirectionalSolver)
    solver.SetMapper(Pool(2).map)
    # solver.SetGenerationMonitor(stepmon)
    solver.SetStrictRanges(min=[0] * FAC_NUM, max=[1] * FAC_NUM)
    solver.SetConstraints(constraint)
    solver.SetPenalty(penalty)
    solver.Solve(helper_cost, termination=COG(1e-07, 15), disp=False)

    weight = solver.Solution()
    func_v = _expect_value(weight, rf, gamma, fac_array)

    values = np.append(arr=weight, values=[func_v, seed, nbins])
    idx = FAC_NAME + ['func_v', 'seed', 'nbins']
    result = pd.Series(data=values, index=idx)
    return result


@click.command()
@click.option('--seed', type=int)
@click.option('--nbins', type=int)
@click.option('--gamma', type=int, default=7)
def main(seed, nbins, gamma):
    equation = equation_str(mr=0.2)
    pf = generate_penalty(generate_conditions(equation), k=1e20)
    cf = generate_constraint(generate_solvers(simplify(equation)))

    def grouped_apply_opt_part(part_df):
        weights_applyed: pd.DataFrame = part_df.groupby('date').apply(
            opti_latti_pow,
            # meta=meta_dict,
            nbins=nbins,
            gamma=gamma,
            constraint=cf,
            penalty=pf,
            seed=seed)
        weights_applyed = weights_applyed.astype({
            'seed': 'i8',
            'nbins': 'i8'
        })
        return weights_applyed

    merged_dd: dd.DataFrame = dd.from_pandas(MERGED_DF, npartitions=2)
    meta_dict = {fac_name: float for fac_name in FAC_NAME + ['func_v']}
    meta_dict.update({'seed': 'i8', 'nbins': 'i8'})

    dd_weights_set = merged_dd.map_partitions(func=grouped_apply_opt_part,
                                              meta=meta_dict)
    with ProgressBar():
        best_weights = dd_weights_set.compute()

    # best_weight: pd.DataFrame = dd_applyed.compute()
    best_weights.to_pickle(
        "data/interim/best_weight_s{}_nb{}_ga{}.pickle".format(
            seed, nbins, gamma))


# %%
if __name__ == "__main__":
    random_num = read_simul_random_df()
    DATE_LIST: list = random_num.index.unique().tolist()
    FAC_NAME: list = random_num.columns.to_list()
    FAC_NUM = random_num.shape[1]

    rf_df = read_rf_df()
    MERGED_DF = join_rf_df(random_num_df=random_num,
                           rf_df=rf_df,
                           rf_type='week')
    if __debug__:
        MERGED_DF = MERGED_DF.loc[DATE_LIST[0:8], :]
    main()
