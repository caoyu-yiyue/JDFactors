import pandas as pd


def read_garch_stand_residual():
    """
    读取由r 语言生成的garch standarized residuals
    """
    garch_resids: pd.DataFrame = pd.read_feather(
        'data/interim/garch_residuals.feather')
    garch_resids.set_index('date', inplace=True)
    return garch_resids


def read_garch_cop_randoms(cop_type: str):
    """
    读取由r 语言生成的garch + copula 模型生成的随机数
    """
    if cop_type == 't':
        rands = pd.read_feather(
            'data/interim/garch_tcop_stMargin_randoms.feather')
    elif cop_type == 'norm':
        rands = pd.read_feather(
            'data/interim/garch_normcop_stMargin_randoms.feather')

    return rands
