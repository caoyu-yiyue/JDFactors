# %%
import numpy as np
import pandas as pd
# from mystic import termination
from mystic.penalty import linear_equality, linear_inequality
from mystic.pools import SerialPool as Pool
from mystic.search import Searcher
from mystic.solvers import BuckshotSolver, PowellDirectionalSolver, diffev2
from mystic.strategy import Rand1Bin
from mystic.symbolic import (generate_conditions, generate_constraint,
                             generate_penalty, generate_solvers, simplify)
from mystic.termination import ChangeOverGeneration as COG

# from mystic.pools import SerialPool as Pool

# from mystic.tools import random_seed

# %%
random_num: pd.DataFrame = pd.read_feather('data/interim/random_num.feather')
random_num['date'] = pd.to_datetime(random_num['date'], format='%Y-%m-%d')
random_num.set_index('date', inplace=True)
grouped = random_num.groupby("date")
first_group: pd.DataFrame = grouped.get_group('2003-01-03')
first_array = first_group.to_numpy()

# %%
rf_df: pd.DataFrame = pd.read_csv('data/raw/csvFiles/TRD_Nrrate.csv',
                                  sep='\t',
                                  parse_dates=['Clsdt'],
                                  index_col='Clsdt',
                                  encoding='utf-16')

# %%
DATE_LIST: list = random_num.index.unique().tolist()
FAC_NUM = random_num.shape[1]

# %%
# 通过现在所在的date 找到下一个date，并在rf 数据框中找到rf 数据
current_date = first_group.index[0]
nxt_date_idx = DATE_LIST.index(current_date) + 1
nxt_date = DATE_LIST[nxt_date_idx]
rf = rf_df.loc[nxt_date, 'Nrrwkdt']


# %%
# @jit(nopython=True)
def expect_value(weight, rf, gamma, fac_ndarr):
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
def weight_sum_cons(weight):
    """权重和为1"""
    return np.sum(weight) - 1


def margin_rate_cons(weight, mr):
    return weight[0] - 2 * (np.sum(weight[1:])) - (1 / mr)


@linear_inequality(margin_rate_cons, kwds={'mr': 0.2}, k=1e20)
@linear_equality(weight_sum_cons, k=1e20)
def penalty(weight):
    return 0.0


# %%
def equation_str(mr):
    equations = """
    x0 + x1 + x2 == 1
    x0 + 2*(x1 + x2) - (1. / {}) <= 0
    """.format(mr)

    return equations


equation = equation_str(mr=0.2)
pf = generate_penalty(generate_conditions(equation), k=1e20)
cf = generate_constraint(generate_solvers(simplify(equation)))
bounds = [(0, 1)] * FAC_NUM
# random_seed(101)

# %%
result = diffev2(cost=expect_value,
                 args=(rf, 7, first_array),
                 x0=bounds,
                 bounds=bounds,
                 gtol=200,
                 scale=0.8,
                 constraint=cf,
                 penalty=pf,
                 strategy=Rand1Bin,
                 full_output=True,
                 npop=800)
result

# %%
# diff_solver = DifferentialEvolutionSolver2(dim=FAC_NUM, NP=40)
# diff_solver.SetConstraints(cf)
# diff_solver.SetPenalty(pf)
# diff_solver.SetStrictRanges(min=[0] * FAC_NUM, max=[1] * FAC_NUM)
# diff_solver.SetInitialPoints(x0=)
# term = termination.CollapseAt(generations=100)
# diff_solver.SetTermination(termination=term)
# diff_solver.SetMapper(Pool().map)
# diff_solver.Solve(cost=expect_value,
#                   ExtraArgs=(rf, 7, first_array),
#                   ScalingFactor=0.7, disp=True)
# diff_solver.Solution()

# %%
stop = COG()
_map = Pool().map


def helper(weight):
    return expect_value(weight, rf=rf, gamma=7, fac_ndarr=first_array)


seacher = Searcher(npts=10,
                   retry=1,
                   tol=8,
                   memtol=1,
                   map=_map,
                   sprayer=BuckshotSolver,
                   seeker=PowellDirectionalSolver)
seacher.Search(model=helper,
               bounds=bounds,
               stop=stop,
               constraints=cf,
               penalty=pf)
