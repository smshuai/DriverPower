import pickle

# Default parameter for XGBoost used in DriverPower
param = {'max_depth':8,
         'eta':0.05,
         'subsample':0.6,
         'nthread':15,
         'objective':'count:poisson',
         'verbose':0,
         'max_delta_step':1.2,
         'eval_metric':'poisson-nloglik'}

# Dump to pickle
with open('xgb.param.pkl', 'wb') as f:
    pickle.dump(param, f)
