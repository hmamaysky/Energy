#!/user/kh3191/.conda/envs/nlp/bin/python
import pandas as pd
import numpy as np
seed = 12345
np.random.seed(seed)
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm

from OOSfuncs import data_set, get_end_of_week, get_test_row_range, PCA_augment, ind_var_list
from tqdm import tqdm
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import Lasso

from datetime import datetime
from dateutil.relativedelta import relativedelta
from glob import glob

import torch

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--wk', type=int, default=8)
    parser.add_argument('--window', type=int, default=5)
    parser.add_argument('--frequency', type=int, default=1)
    parser.add_argument('--cvs', type=int, default=3)
    parser.add_argument('--rolling', type=bool, default=True)
    parser.add_argument('--select_significant', type=bool, default=True)
    parser.add_argument('--top_R2', type=int, default=None)
    opt = parser.parse_args()
    return opt

def standardize(df):
#     if window:
#         rolling_mean = df.rolling(window).mean()
#         rolling_std = df.rolling(window).std()
#         df_scaled = (df - rolling_mean)/ rolling_std
#         return df_scaled.dropna()
#     else:
    return (df - df.mean())/ df.std()

def get_last_month_yyyymm(timestamp):
    # Check if the given date is the last day of the month
    last_day_of_month = timestamp.replace(day=1) + relativedelta(months=1) - relativedelta(days=1)
    if timestamp.day == last_day_of_month.day:
        return timestamp.strftime("%Y%m")
    else:
        # Subtract one month
        last_month = timestamp - relativedelta(months=1)
        return last_month.strftime("%Y%m")

    
def get_train_test_split(d_var, forecast_start, wk, window, rolling, 
                 concat_path='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info',
                ind_vars = ['entropy', 'artcount', 'sCo', 'fCo',
           'sGom', 'fGom', 'sEnv', 'fEnv', 'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc',
           'fRpc', 'sEp', 'fEp', 'PCAsent', 'PCAall', 'PCAfreq']):
    
    data = data_set(d_var)

    ### 1. Get update and forecast window
    date_update_range, date_test_range, date_pca_range=get_test_row_range(data['date'], forecast_start, wk=wk, update_window=window)

    ### 2. Shift x to match the lag 
    lag_vars = ind_var_list(d_var, weeks=wk)
    # trend and WIPIyoy will not lag
    lag_vars.remove('trend')
    lag_vars.remove('WIPI_{}wk'.format(wk))
    data_x = data.copy()
    data_x.loc[:,lag_vars] = data_x.loc[:,lag_vars].shift(wk)

    ### 3. Prepare LHS data for training and testing
    data_ytrain = data[date_update_range]
    data_ytest = data[date_test_range]

    ### 4. Prepare RHS data for training and testing   
    # 4.1 Add PCA series
    data_x_pca = data_x[date_pca_range]

    if not rolling:
        data_x_pca = PCA_augment(data_x_pca)

        # 4.2 Select the train set (first few rows) and the test set (the last row)
        data_xtrain = data_x_pca.iloc[:np.sum(date_update_range),:]

    else:
        timestamp = data_x_pca['date'].max()
        YYYYMM_end = get_last_month_yyyymm(timestamp)

        agg_path_list = glob(f'{concat_path}/*{YYYYMM_end}_NYtime_daily_level_measures_C_2023.csv')
        assert len(agg_path_list) == 1

        df = pd.read_csv(agg_path_list[0]).drop(columns=['Unclassified', 'Unclassified-Sentiment'])
        df['date'] = pd.to_datetime(df['date'].astype(str))
        df.columns = df.columns.str.replace('Topic-Sentiment ', 's')
        df.columns = df.columns.str.replace('Topic ', 'f')
        df.rename(columns={'article count': 'artcount'}, inplace=True)
        df.set_index('date', inplace=True)
        df = df.resample(f'W-{get_end_of_week(d_var).upper()}').mean().rolling(4).mean()
        df = df.shift(wk)

        data_x_pca = pd.merge(data_x_pca[['date']], df, on='date', how='left')
        data_x_pca.set_index('date', inplace=True)
        data_x_pca = PCA_augment(standardize(data_x_pca).dropna(), topic_labels=False)
        data_x_pca = standardize(data_x_pca)

        # 4.2 Select the train set (first few rows) and the test set (the last row)
        data_ytrain = pd.merge(data_ytrain, data_x_pca, on='date')
        data_xtrain = pd.merge(data[date_update_range]['date'], data_x_pca, on='date')


    assert len(data_xtrain) == len(data_ytrain)

    data_xtest = data_x_pca.iloc[[-1],:]
    assert len(data_xtest) == len(data_ytest)

    ### 5. Lasso Update and Prediction Process
    # 5.1 First, use grid search to choose the best Penalty 
    #     Then, run lasso using that panalty regression and update, 
    #     here, sm.add_constant adds a constant to RHS
    if not rolling:
        x_train = data_xtrain.loc[:,ind_vars]
    else:
        x_train = data_xtrain.drop(columns=['date'])

    ## Choose the correct LHS var according to forecasting duration
    y_train = data_ytrain[f'{d_var}_t{wk}']

    scaler = {'ymean': y_train.mean(), 'ystd': y_train.std()}
    train_xy = pd.concat([x_train,y_train],axis=1).dropna()
    # scaling is necessary in LASSO
    train_xy = standardize(train_xy)
    y_train = train_xy.iloc[:,-1]
    x_train = train_xy.iloc[:,0:-1]
    
    return x_train, data_xtest, y_train, data_ytest, YYYYMM_end, scaler
    
if __name__ == '__main__':
    
    opt = parse_option()
    print(opt)

    for d_var in ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']:

        print(d_var)

        time_col = data_set(d_var)['date'] 
        time_lower = time_col[0]
        # time_lower= pd.Timestamp('2014-05-03')
        test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*opt.window*52)+'days')][::opt.frequency]

        best_lambda_dic = {}
        significant_ind_vars_dic = {}
        
        for ii, forecast_start in enumerate(tqdm(test_week_list)): #[876:]

            try:
                x_train, data_xtest, y_train, data_ytest, YYYYMM_end, scaler = get_train_test_split(d_var, 
                                                                                            forecast_start, 
                                                                                            opt.wk, 
                                                                                            opt.window, 
                                                                                            opt.rolling)
                
                ## alpha = 0: Manual Cross-Validation for OLS
                kf = KFold(n_splits=opt.cvs)

                if opt.select_significant:
                    avg_r2_score = {}
                    for d_var_one_predictor in x_train.columns:
                        x_train_one_predictor = x_train[d_var_one_predictor]

                        r2_scores = []
                        for train_index, test_index in kf.split(x_train_one_predictor):
                            X_cv_train, X_cv_test = x_train_one_predictor.iloc[train_index], x_train_one_predictor.iloc[test_index]
                            y_cv_train, y_cv_test = y_train.iloc[train_index], y_train.iloc[test_index]
                            model = lm.OLS(y_cv_train, X_cv_train, missing='drop').fit()
                            y_cv_pred = model.predict(X_cv_test)#sm.add_constant(X_cv_test))
                            r2_scores.append(r2_score(y_cv_test, y_cv_pred))

                        avg_r2_score[d_var_one_predictor] = np.mean(r2_scores)
                    significant_ind_vars = [key for key, value in sorted(avg_r2_score.items(), key=lambda item: item[1], reverse=True)[:opt.top_R2]]
                    significant_ind_vars_dic[str(forecast_start)[:10]] = significant_ind_vars
                    
                    x_train = x_train[significant_ind_vars]
                    data_xtest = data_xtest[significant_ind_vars]
                    

#                 mse_scores = []
#                 for train_index, test_index in kf.split(x_train):
#                     X_train = x_train#sm.add_constant(x_train)
#                     X_cv_train, X_cv_test = X_train.iloc[train_index], X_train.iloc[test_index]
#                     y_cv_train, y_cv_test = y_train.iloc[train_index], y_train.iloc[test_index]
#                     model = lm.OLS(y_cv_train, X_cv_train, missing='drop').fit()
#                     y_cv_pred = model.predict(X_cv_test)#sm.add_constant(X_cv_test))
#                     mse_scores.append(mean_squared_error(y_cv_test, y_cv_pred))
#                 neg_avg_mse = -np.mean(mse_scores)
                neg_avg_mse = -99999999

                ## Set up Lasso instance and grid search for penalty coefficient
                pre_model = Lasso(random_state=seed, fit_intercept=False)
                param_grid = [{'alpha': np.exp(np.linspace(-7, 0, 40))}]
                grid_search = GridSearchCV(pre_model, param_grid, cv=opt.cvs, scoring='neg_mean_squared_error', n_jobs=-1)
                grid_search.fit(x_train, y_train)

                if neg_avg_mse > grid_search.best_score_:
                    best_lambda = 0
                else:
                    best_lambda = grid_search.best_params_['alpha']

                    
                if not opt.rolling:
                    best_lambda_dic[str(forecast_start)[:10]] = best_lambda
                else:
                    best_lambda_dic[str(forecast_start)[:10]] = [YYYYMM_end, best_lambda, scaler]

            except AssertionError:
                break


        if not opt.rolling:
            torch.save(best_lambda_dic, f'res_all_variables/{opt.cvs}fold/forward_best_lambda_{d_var}.pt')
            
        else:
            torch.save(best_lambda_dic, f'res_all_variables/{opt.cvs}fold/forward_rolling_best_lambda_{d_var}.pt')
            if opt.select_significant:
                torch.save(significant_ind_vars_dic, f'res_all_variables/{opt.cvs}fold/forward_rolling_selected_{d_var}.pt')
