# %%
# import dask.dataframe as dd
# from dask.diagnostics import ProgressBar
import numpy as np
import pandas as pd
import click
from tqdm import tqdm
import os
from mystic.solvers import (LatticeSolver, PowellDirectionalSolver,
                            DifferentialEvolutionSolver2)
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
    # r_{f,t+1} 与r_{t+1} 是同一时间刻度上，所以只需要简单对齐合并即可
    shifted_rf: pd.Series = rf_df[rf_type]
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


def equation_str(mr, sum_1):
    if sum_1:
        equations = """
        x0 + x1 + x2 == 1
        x0 + 2*(x1 + x2) - (1. / {}) <= 0
        """.format(mr)
    else:
        equations = """
        x0 + 2*(x1 + x2) - (1. / {}) <= 0
        """.format(mr)

    return equations


def opti_fun(df: pd.DataFrame, gamma: int, max_r: float, seed: int, nbins: int,
             method: str, constraint, penalty, FAC_NAME, FAC_NUM):
    # print(df.name)
    core_fac = df.drop(columns='rf')
    fac_array = core_fac.to_numpy()
    # fac_names = core_fac.columns
    rf = df['rf'][0]
    random_seed(s=seed)

    def helper_cost(weight):
        return _expect_value(weight, rf=rf, gamma=7, fac_ndarr=fac_array)

    if method == 'Lattice':
        solver = LatticeSolver(dim=FAC_NUM, nbins=nbins)
        solver.SetNestedSolver(PowellDirectionalSolver)
        solver.SetMapper(Pool(2).map)
    elif method == 'DE':
        solver = DifferentialEvolutionSolver2(dim=FAC_NUM, NP=nbins)
        solver.SetRandomInitialPoints(min=[0] * FAC_NUM, max=[max_r] * FAC_NUM)
    # solver.SetGenerationMonitor(stepmon)
    solver.SetStrictRanges(min=[0] * FAC_NUM, max=[max_r] * FAC_NUM)
    solver.SetConstraints(constraint)
    solver.SetPenalty(penalty)

    if method == 'Lattice':
        solver.Solve(helper_cost, termination=COG(1e-07, 15), disp=False)
    elif method == 'DE':
        solver.Solve(helper_cost,
                     termination=COG(1e-07, 60),
                     disp=False,
                     ScalingFactor=0.7)

    weight = solver.Solution()
    func_v = solver.bestEnergy

    values = np.append(arr=weight, values=[func_v])
    idx = FAC_NAME + ['func_v']
    result = pd.Series(data=values, index=idx)
    return result


@click.command()
@click.option('--gamma', type=int, default=7)
@click.option('--mr', type=float)
@click.option('--max_r', type=str, default='None')
@click.option('--sum_1/--no_sum_1', default=True)
@click.option('--seed', type=int)
@click.option('--nbins', type=int)
@click.option('--method', type=click.Choice(['Lattice', 'DE']))
@click.option('--half', type=click.Choice(['first', 'second']))
@click.argument('output_file', type=click.Path(writable=True, dir_okay=True))
def main(gamma, mr, max_r, sum_1, seed, nbins, method, half, output_file):
    day_len = len(DATE_LIST)
    split_point = day_len // 2
    if half == 'first':
        half_dates = DATE_LIST[:split_point]
    elif half == 'second':
        half_dates = DATE_LIST[split_point:]

    if __debug__:
        use_df = RANDOM_NUM.loc[half_dates[0:8], ]
    else:
        use_df = RANDOM_NUM.loc[half_dates, :]

    rf_df = read_rf_df()
    merged_df = join_rf_df(random_num_df=use_df, rf_df=rf_df, rf_type='week')

    equation = equation_str(mr=mr, sum_1=sum_1)

    if max_r == 'None':
        max_r = None
    else:
        max_r = float(max_r)

    pf = generate_penalty(generate_conditions(equation), k=1e20)
    cf = generate_constraint(generate_solvers(simplify(equation)))

    tqdm.pandas()
    weights_applyed: pd.DataFrame = merged_df.groupby('date').progress_apply(
        opti_fun,
        # meta=meta_dict,
        gamma=gamma,
        max_r=max_r,
        seed=seed,
        nbins=nbins,
        method=method,
        constraint=cf,
        penalty=pf,
        FAC_NAME=FAC_NAME,
        FAC_NUM=FAC_NUM)
    # 为每一行追加seed、nbins、method 信息。先开一个meta_df，然后与所需数据合并
    meta_df: pd.DataFrame = pd.DataFrame(data={
        'seed': seed,
        'nbins': nbins,
        'method': method
    },
                                         index=weights_applyed.index)
    weights_with_meta: pd.DataFrame = pd.concat([weights_applyed, meta_df],
                                                axis=1)

    # 将seed 和nbins 转化成int8 类型
    weights_with_meta = weights_with_meta.astype({
        'seed': 'int8',
        'nbins': 'int8'
    })
    # 使用r_{f, t+1} 和r_{t+1} 计算的是前一个时刻的权重w_{t}，所以结果需要向前提一个时刻
    # weights_adj_date = weights_applyed.shift(-1).dropna()

    # merged_dd: dd.DataFrame = dd.from_pandas(merged_df, npartitions=2)
    # meta_dict = {fac_name: float for fac_name in FAC_NAME + ['func_v']}
    # meta_dict.update({'seed': 'i8', 'nbins': 'i8'})

    # dd_weights_set = merged_dd.map_partitions(func=grouped_apply_opt_part,
    #                                           meta=meta_dict)
    # start = time.time()
    # with ProgressBar():
    #     best_weights = dd_weights_set.compute()
    # stop = time.time()
    # print(stop - start)

    # best_weight: pd.DataFrame = dd_applyed.compute()

    # 如果保存文件的路径不存在，则创建一个
    save_dir = os.path.split(output_file)[0]
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    weights_with_meta.to_pickle(output_file)


# %%
if __name__ == "__main__":
    RANDOM_NUM = read_simul_random_df()
    DATE_LIST: list = RANDOM_NUM.index.unique().tolist()
    FAC_NAME: list = RANDOM_NUM.columns.to_list()
    FAC_NUM = RANDOM_NUM.shape[1]

    main()
