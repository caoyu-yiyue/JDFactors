import pandas as pd


def read_facs(fpath='data/raw/csvFiles/week_readed.csv'):
    """
    读取准备好的数据框
    Return:
    -------
    pd.DataFrame
        整理好的数据框
    """
    facs_ret = pd.read_csv(fpath, index_col='trdWeek', parse_dates=['trdWeek'])
    facs_ret.rename_axis(index={'trdWeek': 'date'}, inplace=True)
    return facs_ret.loc['2001-01-01':]


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


def describe_df(df_for_describe: pd.DataFrame):
    """
    传入一个 DataFrame 或 Series，在描述性统计的基础上加一个偏度和峰度
    Return:
    ------
    pd.DataFrame
    """
    des_df = df_for_describe.describe()
    skew_v = df_for_describe.skew()
    kurt_v = df_for_describe.kurtosis()
    if isinstance(des_df, pd.DataFrame):
        des_df = des_df.append(skew_v.rename('skew'))
        des_df = des_df.append(kurt_v.rename('kurtosis'))
    elif isinstance(des_df, pd.Series):
        skew_kurt: pd.Series = pd.Series(data=[skew_v, kurt_v],
                                         index=['skew', 'kurt'])
        des_df = des_df.append(skew_kurt)
    return des_df
