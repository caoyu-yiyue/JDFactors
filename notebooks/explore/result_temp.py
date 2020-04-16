# %%
from src.data import preparing_data as preda
from src.process import cum_return as pcr
import pandas as pd
from IPython.core.interactiveshell import InteractiveShell
import matplotlib.pyplot as plt
InteractiveShell.ast_node_interactivity = 'all'

# %%
facs: pd.DataFrame = preda.read_facs()
rf_week = preda.read_rf_df()['week']

# %%
param_dicts = [{
    'gamma': 7,
    'mr': 0.2,
    'max_r': 1,
    'sum1': True
}, {
    'gamma': 7,
    'mr': 0.5,
    'max_r': 1,
    'sum1': True
}, {
    'gamma': 3,
    'mr': 0.5,
    'max_r': 1,
    'sum1': True
}]

# %%
for param_dict in param_dicts:
    print(param_dict)
    ext_ret = pcr.ex_ret_by_weight(gamma=param_dict['gamma'],
                                   mr=param_dict['mr'],
                                   max_r=param_dict['max_r'],
                                   sum1=param_dict['sum1'],
                                   fac_returns=facs)
    full_ret = ext_ret.add(rf_week).dropna()
    preda.describe_df(full_ret)
    full_ret.hist()
    plt.show()
    full_ret.plot()
    plt.show()
