"""一些描述性统计，观察数据结构"""

# %%
import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.tsa.api as smt
import matplotlib.pyplot as plt
import itertools
# from arch.univariate import arch_model
from IPython.core.interactiveshell import InteractiveShell

# %%
# 设置全部显示
InteractiveShell.ast_node_interactivity = 'all'

# %%
# read data and prepare index
week_3fac_raw: pd.DataFrame = pd.read_csv(
    'data/raw/csvFiles/ThrfacWeek.csv',
    sep='\t',
    header=0,
    names=['marketTypeID', 'trdWeek', 'premium', 'smb', 'hml'],
    usecols=['trdWeek', 'premium', 'smb', 'hml'],
    encoding='utf-16',
    index_col='trdWeek',
    parse_dates=True)
week_3fac_raw.index = pd.to_datetime(week_3fac_raw.index + '-5',
                                     format='%Y-%W-%w')
week_3fac_raw.sort_index(inplace=True)

# %%
# 选择所需的时间区间
week_3fac: pd.DataFrame = week_3fac_raw.loc['2000-01-01':]

# %%
# 描述性统计
week_3fac.describe()
week_3fac.skew()
week_3fac.kurtosis()
week_3fac.corr()
week_3fac.hist()
week_3fac.plot(subplots=True)

# 因子累积收益率图
((week_3fac + 1).cumprod() - 1).plot()

# %%
# qqplot
for col in week_3fac.columns:
    sm.qqplot(week_3fac[col], line='s', ylabel=col + ' sample quantiles')

# %%
# 散点图矩阵
pd.plotting.scatter_matrix(week_3fac)


# %%
def threshhold_corr(df: pd.DataFrame):
    """
    计算threshhold correlation 的函数

    Parameters:
    -----------
    df: pd.DataFrame
        传入一个两列的DataFrame，返回使用该两列计算的threshhold corr

    Return:
    -------
        tuple of dict
        第一个元素为使用极端小值计算的corr，第二个元素为使用极端大值计算的corr；
        每个元素为一个dict，key 为使用的threshhold（概率值），value 为对应的分位数。
    """
    bottom_left = {}
    up_right = {}

    for u in np.arange(0, 0.5, step=0.01):
        # 两列都小于其各自的u 分位数时，作为有效下极限数值
        extremes: pd.DataFrame = df[
            (df.iloc[:, 0] < df.iloc[:, 0].quantile(u))
            & (df.iloc[:, 1] < df.iloc[:, 1].quantile(u))]
        # 如果极限数值小于20 行，则直接进入下一个循环
        if extremes.shape[0] < 20:
            continue

        extr_corr = extremes.corr().iloc[0, 1]
        bottom_left[u] = extr_corr

    for u in np.arange(0.5, 1, step=0.01):
        # 两列都大于其各自u 分位数时，作为有效的上极限数值
        extremes: pd.DataFrame = df[
            (df.iloc[:, 0] > df.iloc[:, 0].quantile(u))
            & (df.iloc[:, 1] > df.iloc[:, 1].quantile(u))]
        # 如果极限数值小于20 行，直接结束循环
        if extremes.shape[0] < 20:
            break

        extr_corr = extremes.corr().iloc[0, 1]
        up_right[u] = extr_corr

    return bottom_left, up_right


# %%
# 为每个因子的两两组合计算threshhold corr 并画图
for col_pair in itertools.combinations(week_3fac.columns, 2):
    thresh_corr: tuple = threshhold_corr(week_3fac[list(col_pair)])
    plt.plot(*zip(*thresh_corr[0].items()))
    plt.plot(*zip(*thresh_corr[1].items()))
    plt.title('Threshhold Corr for ' + ' and '.join(col_pair))
    plt.show()

# %% [markdown]
# ## 自相关
# 为每个因子及其绝对值画自相关图
for col in week_3fac.columns:
    smt.graphics.plot_acf(week_3fac[col],
                          lags=50,
                          title='Autocorrelation for ' + col)
    smt.graphics.plot_acf(week_3fac[col].abs(),
                          lags=50,
                          title='Autocorrelation for absolute ' + col)

# %%
# 前三个自相关系数
week_3fac.apply(lambda col: smt.stattools.acf(col, qstat=True)[0])[1:4]
# 对应的p 值
week_3fac.apply(lambda col: smt.stattools.acf(col, qstat=True)[2])[0:3].lt(
    0.05)


# %%
# 找到一个最优的ARIMA
def _get_best_ARIMA(ts_series):
    """
    根据AIC 最小的原则，选择最优的ARIMA 模型。p, q 在0～4，d 为0～1。
    Parameters:
    -----------
        ts_series: pd.Series
        需要计算的单维时间序列

    Returns:
    --------
        tuple:
        三项分别为：最优的aic，最优的order，最优的模型对象

    """
    best_aic = np.inf
    best_order = None
    best_mdl = None

    pq_rng = range(5)
    d_rng = range(2)

    for i in pq_rng:
        for d in d_rng:
            for j in pq_rng:
                try:
                    tmp_mdl = smt.ARIMA(ts_series,
                                        order=(i, d, j)).fit(method='mle',
                                                             trend='nc')
                    tmp_aic = tmp_mdl.aic
                    if tmp_aic < best_aic:
                        best_aic = tmp_aic
                        best_order = (i, d, j)
                        best_mdl = tmp_mdl
                except Exception:
                    continue

    print('Best AIC: {:6.5f} | order: {}'.format(best_aic, best_order))
    return [best_aic, best_order, best_mdl]


arima_3fac = week_3fac.apply(_get_best_ARIMA, result_type='expand')
"""
best_orders
premium    (2, 0, 1)
smb        (4, 0, 4)
hml        (4, 0, 2)
"""

# %%
# resid 序列的acf pcf
for col in arima_3fac.columns:
    smt.graphics.plot_acf(arima_3fac.loc[2, col].resid,
                          title='Residual Autocorr for ' + col)
    smt.graphics.plot_pacf(arima_3fac.loc[2, col].resid,
                           title='Residual Partial Autocorr ' + col)

# %%
# 查看统一模型(2, 0, 1)的相关系数
one_mdl_smb = smt.ARIMA(week_3fac['smb'], order=(2, 0, 1)).fit(method='mle',
                                                               trend='nc')
one_mdl_hml = smt.ARIMA(week_3fac['hml'], order=(2, 0, 1)).fit(method='mle',
                                                               trend='nc')

# %%
smt.graphics.plot_acf(one_mdl_smb.resid,
                      title='Residual Autocorr for smb (2, 0, 1)')
smt.graphics.plot_pacf(one_mdl_smb.resid,
                       title='Residual Partial Autocorr for smb (2, 0, 1)')

# %%
smt.graphics.plot_acf(one_mdl_hml.resid,
                      title='Residual Autocorr for hml (2, 0, 1)')
smt.graphics.plot_pacf(one_mdl_hml.resid,
                       title='Residual Partial Autocorr for hml (2, 0, 1)')

# %%
# 使用garch 模型
# am = arch_model(100 * week_3fac['hml'], mean='AR', lags=6, o=1, dist='skewt')
# res = am.fit()
# smt.graphics.plot_acf(res.resid.dropna() / 100, lags=50)

# smt.graphics.plot_acf(res.resid.dropna().abs() / 100, lags=50)
