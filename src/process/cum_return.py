# %%
import pandas as pd


def read_facs(fpath='data/raw/csvFiles/week_readed.csv'):
    facs_ret = pd.read_csv(fpath, index_col='trdWeek', parse_dates=['trdWeek'])
    facs_ret.rename_axis(index={'trdWeek': 'date'}, inplace=True)
    return facs_ret


def cum_return(gamma: int, max_r: float, sum1: bool,
               fac_returns: pd.DataFrame):
    """
    输入gamma max_r 和sum_1，返回计算所得的累积收益率

    Paramters:
    ----------
    gamma: int
        计算最优化权重所用的gamma 值
    max_r: folat
        计算最优权重时，权重的最大值
    sum1: Bool
        计算最优权重时，权重和是否等于1.
    fac_returns: pd.DataFrame
        因子回报的序列

    Return:
    -------
    folat
        根据权重值与因子收益计算得到的时序累积收益率
    """

    sum1_str = 'T' if sum1 else 'F'
    data_path = ('data/processed/ga{}_max{}_sum1{}/'
                 'combined_result.pickle'.format(gamma, max_r, sum1_str))
    weights_df = pd.read_pickle(data_path)

    fac_names = fac_returns.columns.to_list()
    facs_ret_aligned = fac_returns.reindex(weights_df.index)
    weights_df_aligned: pd.DataFrame = weights_df[fac_names]

    ret_series: pd.Series = weights_df_aligned.mul(
        facs_ret_aligned, axis='index').sum(axis='columns')

    cum_ret = (ret_series + 1).prod() - 1
    return cum_ret
