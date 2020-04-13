# %%
import pandas as pd


def ex_ret_by_weight(gamma: int, mr: float, max_r: float, sum1: bool,
                     fac_returns: pd.DataFrame):
    """
    输入gamma max_r 和sum_1，返回计算所得的累积收益率

    Paramters:
    ----------
    gamma: int
        计算最优化权重所用的gamma 值

    mr: float
        计算最优化权重时使用的mr 值

    max_r: folat
        计算最优权重时，权重的最大值

    sum1: Bool
        计算最优权重时，权重和是否等于1.

    fac_returns: pd.DataFrame
        因子回报的序列

    Return:
    -------
    pd.Series
        根据权重值与因子收益计算得到的没时刻收益率序列
    """

    sum1_str = 'T' if sum1 else 'F'
    data_path = ('data/processed/ga{}_mr{}_max{}_sum1{}/'
                 'combined_result.pickle'.format(gamma, mr, max_r, sum1_str))
    weights_df = pd.read_pickle(data_path)

    fac_names = fac_returns.columns.to_list()
    facs_ret_aligned = fac_returns.reindex(weights_df.index)
    weights_df_aligned: pd.DataFrame = weights_df[fac_names].astype(float)

    ex_ret_series: pd.Series = weights_df_aligned.mul(
        facs_ret_aligned, axis='index').sum(axis='columns')

    return ex_ret_series


def cum_ret(ret_series, rf_series=None):
    """
    输入一列收益率序列，一列可选的rf_series，返回历史收益率

    Parameters:
    -----------
    ret_series: pd.Series
        收益率序列

    rf_series: pd.Series, default None
        可选的无风险收益率序列

    Return:
        pd.Series
        历史累积收益率序列
    """
    if rf_series:
        ret_series = ret_series + rf_series
    cum_ret = (ret_series + 1).cumprod() - 1
    return cum_ret


def annualized_ret(ret_series: pd.Series, feq: str = 'week'):
    """
    计算年化收益率的函数

    Paramters:
    ----------
    ret_series: pd.Series
        收益率序列

    feq: str, one of ('week', 'day', 'month')
        原收益率序列的频率

    Return:
    -------
        float
        年化累积收益率
    """
    if feq == 'week':
        one_year_len = 52
    elif feq == 'day':
        one_year_len = 365
    elif feq == 'month':
        one_year_len = 12

    held_len = len(ret_series)
    annua_ret = ret_series.add(1).prod()**(one_year_len / held_len) - 1
    return annua_ret
