# %%
import pandas as pd
import numpy as np
import itertools
from matplotlib import pyplot as plt
import statsmodels.tsa.api as smt


def _threshhold_corr(df: pd.DataFrame):
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


def threshhold_corr_plot(facs_dfs: dict, colors, binorm_rand=None):
    if binorm_rand is not None:
        binorm_thresh_corr = _threshhold_corr(binorm_rand)
    for col_pair in itertools.combinations(facs_dfs.columns, 2):
        thresh_corr: tuple = _threshhold_corr(facs_dfs[list(col_pair)])
        plt.plot(*zip(*thresh_corr[0].items()),
                 color=colors[0],
                 label='factor_ret')
        plt.plot(*zip(*thresh_corr[1].items()), color=colors[0])
        if binorm_rand is not None:
            plt.plot(*zip(*binorm_thresh_corr[0].items()),
                     color=colors[1],
                     label='norm')
            plt.plot(*zip(*binorm_thresh_corr[1].items()), color=colors[1])
        plt.title('Threshhold Corr for ' + ' and '.join(col_pair))
        plt.legend()
        plt.show()


def autocorr_func_plot(fac_df: pd.DataFrame):
    fig, axs = plt.subplots(nrows=fac_df.shape[1], ncols=2, figsize=[10, 8])
    plt.subplots_adjust(hspace=0.5)
    for i, col in enumerate(fac_df.columns):
        smt.graphics.plot_acf(fac_df[col],
                              ax=axs[i, 0],
                              lags=50,
                              title='Autocorrelation for ' + col)
        smt.graphics.plot_acf(fac_df[col].abs(),
                              lags=50,
                              ax=axs[i, 1],
                              title='Autocorrelation for absolute ' + col)


# %%
# if __name__ == "__main__":
#     facs = preda.read_facs()
#     multinorm_ran = np.random.multivariate_normal(mean=[0, 0],
#                                                   cov=[[1, 0], [0, 1]],
#                                                   size=10000)
#     multinorm_ran_df: pd.DataFrame = pd.DataFrame(multinorm_ran)
#     threshhold_corr_plot(facs,
#                          colors=['red', 'green'],
#                          binorm_rand=multinorm_ran_df)
