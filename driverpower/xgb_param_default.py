import pickle

# Default parameter for XGBoost used in DriverPower
param = {'max_depth': 8,
         'eta': 0.05,
         'subsample': 0.6,
         'nthread': 15,
         'objective': 'count:poisson',
         'max_delta_step': 1.2,
         'eval_metric': 'poisson-nloglik',
         'silent': 1,
         'verbose_eval': 100,  # print evalutation every 100 rounds
         'early_stopping_rounds': 5,
         'num_boost_round': 5000  # max number of rounds
        }

# Dump to pickle
with open('xgb_param.pkl', 'wb') as f:
    pickle.dump(param, f)