#!/user/kh3191/.conda/envs/nlp/bin/python
from OOSfuncs import data_set, get_test_row_range, ind_var_list, PCA_augment, standardize, get_end_of_week, get_end_of_week_for_rolling_average
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold, GridSearchCV
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import Lasso
seed = 12345
import os
import torch
from tqdm import tqdm

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('--wk', type=int, default=8,
                        help='Number of weeks for the shift operation. Must be an integer.')

    parser.add_argument('--window', type=int, default=5,
                        help='Size of the window for rolling operations. Must be an integer.')

    parser.add_argument('--frequency', type=int, default=1,
                        help='Frequency for data processing. Must be an integer.')

    parser.add_argument('--cvs', type=int, default=3,
                        help='Number of cross-validation splits. Must be an integer.')

    parser.add_argument('--rolling', type=bool, default=True,
                        help='Flag to indicate if rolling operation should be applied. True or False.')

    parser.add_argument('--select_significant', type=bool, default=True,
                        help='Flag to select significant features. True or False.')
    
    parser.add_argument('--num_nontext_var', type=int, default=0,
                        help='Number of top R-squared values to consider for non-text variables. Must be an integer.')
    
    parser.add_argument('--num_text_var', type=int, default=2,
                        help='Number of top R-squared values to consider for text varialbes. Must be an integer.')

    parser.add_argument('--normalize', type=bool, default=True,
                        help='Flag to normalize data. True or False.')

    parser.add_argument('--model_name', type=str, default='Forward_Lasso',
                        help='Name of the model to be used. Must be a string.')

    parser.add_argument('--d_var', type=str, default='FutRet', 
                        choices=['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd'],
                        help='Dependent variable to be used in the model. Choose from predefined options.')
    
    opt = parser.parse_args()
    return opt

from dateutil.relativedelta import relativedelta
# def get_last_month_yyyymm(timestamp):
#     # Check if the given date is the last day of the month
#     last_day_of_month = timestamp.replace(day=1) + relativedelta(months=1) - relativedelta(days=1)
#     if timestamp.day == last_day_of_month.day:
#         return timestamp.strftime("%Y%m")
#     else:
#         # Subtract one month
#         last_month = timestamp - relativedelta(months=1)
#         return last_month.strftime("%Y%m")
def get_last_month_yyyymm(timestamp):
    # changed on Feb 27 2024, because end_of_week=Fri means the observations are on Thursday, so should not average on rolling basis on Fridays
    last_month = timestamp - relativedelta(months=1)
    return last_month.strftime("%Y%m")

def adjust_to_friday(date):
    if date.weekday() == 4:
        return date
    else:
        return date + pd.Timedelta(days=(4-date.weekday()) % 7)
    
def drop_na_xy(x_train, y_train):
    assert (x_train.index == y_train.index).min()
    train_xy = pd.concat([x_train, y_train], axis=1).dropna()
    # split the concatenated dataframe back into features and target variable for training
    x_train = train_xy.iloc[:,:-1]
    y_train = train_xy.iloc[:,-1]
    return x_train, y_train

def fill_na_train_test(x_train, x_test):
    assert (x_train.columns == x_test.columns).min()
    x_train_test = pd.concat([x_train, x_test], axis=0).ffill()
    # split the concatenated dataframe back into features and target variable for training
    x_train = x_train_test.iloc[:-1]
    x_test = x_train_test.iloc[[-1]]
    return x_train, x_test

def get_daily_rolling_topic_series(opt, YYYYMM_topic_end, 
      concat_path='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/rolling_combined_info'):
    
    from glob import glob
    agg_path_list = glob(f'{concat_path}/*{YYYYMM_topic_end}_NYtime_daily_level_measures_C_2023.csv')
    assert len(agg_path_list) == 1

    df = pd.read_csv(agg_path_list[0]).drop(columns=['Unclassified', 'Unclassified-Sentiment'])
    df['date'] = pd.to_datetime(df['date'].astype(str))
    df.columns = df.columns.str.replace('Topic-Sentiment ', 's')
    df.columns = df.columns.str.replace('Topic ', 'f')
    df.rename(columns={'article count': 'artcount'}, inplace=True)
    df.set_index('date', inplace=True)
    df.index.name = 'true_date'
    return df

def get_train_test_split_for_text(opt, end_of_train):
    """
    This function splits a dataset into training and testing sets based on the provided end date for training. 
    It is designed to work with time-series data, particularly for scenarios where train and test sets 
    are temporally separated.

    Inputs:
    - opt: A dictionary or similar structure containing all hyperparameters. 
           This includes the number of folds for cross-validation, flags for enabling/disabling 
           normalization, handling rolling topics, etc.
    - end_of_train: A date string (YYYY-MM-DD) that specifies the last date to be included in the training set.

    Outputs:
    - data_x_train: A pandas DataFrame containing the features for the training set.
    - data_y_train: A pandas Series containing the target variable for the training set.
    - data_x_test: A pandas DataFrame containing the features for the testing set.
    - data_y_test: A pandas Series containing the target variable for the testing set.

    Example Usage:
    Suppose 'end_of_train' is set to '2003-04-25'. This would result in the following splits:
    - x_train: Data from the period [1998-06-26 to 2003-04-25] with columns like 'artcount', 'entropy', 'f*', 's*', 'PCA*'.
    - y_train: Target variable data from the same period [1998-06-26 to 2003-04-25] with the column 'd_var_t8_fwd'.
    - x_test: Data for the specific date [2003-06-20] with columns similar to the training set.
    - y_test: Target variable data for the same date [2003-06-20].
    
    Note: The asterisks (*) in column names like 'f*' and 's*' represent wildcard characters, 
    indicating multiple columns that follow a naming pattern.
    """
    # read in full dataset according to the specified dependent variable (d_var)
    data = data_set(opt.d_var)
    # create forward-looking dependent variables by shifting the specified weeks (wk) forward
    data[f'{opt.d_var}_t{opt.wk}_fwd'] = data[f'{opt.d_var}_t{opt.wk}'].shift(-opt.wk)
    # calculate ranges for updating the dataset, testing, and PCA based on the training end date
    date_update_range, date_test_range, date_pca_range = get_test_row_range(data['true_date'], end_of_train, wk=opt.wk, update_window=opt.window)
    # determine the start and end dates for the training data
    X_train_start_date = data['true_date'][date_update_range].min()
    X_train_end_date = data['true_date'][date_update_range].max()
    assert X_train_end_date == end_of_train

    # get the last month in the format YYYYMM for the beginning of the test range
    YYYYMM_topic_end = get_last_month_yyyymm(data['true_date'][date_test_range].iloc[0])
    # print(YYYYMM_topic_end)
    
    # process data for rolling topics if the 'rolling' option is enabled
    if opt.rolling:
        # get the series for rolling topics up to the specified end month
        df_daily_topic = get_daily_rolling_topic_series(opt, YYYYMM_topic_end)
        df = df_daily_topic.resample(f'W-{get_end_of_week_for_rolling_average(opt.d_var)}').mean().rolling(4).mean()
        # added on Mar 5: label the Thursdays using Fridays
        if get_end_of_week(opt.d_var) == 'Fri':
            df.index += pd.Timedelta(days=1)
        # join additional rows from the PCA range to the rolling topic data, handling missing keys if necessary
        df = df.join(data.loc[date_pca_range, ['true_date']].set_index('true_date'), on='true_date', how='left')
        # augment the data with PCA components and fill forward missing values
        pca = PCA_augment(df.dropna(), topic_labels=False)[['PCAsent', 'PCAfreq', 'PCAall']]
        data_x = df.join(pca, how='left').ffill()
        
        # added on Mar 5
        df = df_daily_topic.resample(f"W-{df_daily_topic.index[-1].strftime('%A')[:3]}").mean().rolling(4).mean()
        # added on Mar 5: label the Thursdays using Fridays
        df.index = df.index.map(adjust_to_friday)
        # join additional rows from the PCA range to the rolling topic data, handling missing keys if necessary
        df = df.join(data.loc[date_pca_range, ['true_date']].set_index('true_date'), on='true_date', how='left')
        # augment the data with PCA components and fill forward missing values
        pca = PCA_augment(df.dropna(), topic_labels=False)[['PCAsent', 'PCAfreq', 'PCAall']]
        data_x_most_recent_topic_for_testing = df.join(pca, how='left').ffill()
        
        
    # process non-rolling data by selecting specific columns for the PCA range
    else:
        data_x = data.loc[date_pca_range, ['true_date', 'artcount', 'entropy', 'sCo', 'fCo', 'sGom', 'fGom', 'sEnv', 'fEnv',
                  'sEpg', 'fEpg', 'sBbl', 'fBbl', 'sRpc', 'fRpc', 'sEp', 'fEp'] + ['PCAsent', 'PCAfreq', 'PCAall']]\
                    .set_index('true_date')
    
    # extract training data based on the start and end dates, sorting by date
    data_x_train = data_x[(data_x.index >= X_train_start_date) & (data_x.index <= X_train_end_date)]\
                    .sort_values('true_date')
    # get the training target variable corresponding to the training feature dates
    if get_end_of_week_for_rolling_average(opt.d_var) == 'Tue':
        data_y_train = data.set_index('date_Tue').loc[data_x_train.index, f'{opt.d_var}_t{opt.wk}_fwd']
    elif get_end_of_week_for_rolling_average(opt.d_var) == 'Thu':
        data_y_train = data.set_index('true_date').loc[data_x_train.index, f'{opt.d_var}_t{opt.wk}_fwd']
    else:
        raise ValueError

    # extract testing data (the last row of data_x)
    data_x_test = data_x_most_recent_topic_for_testing.iloc[[-1]]
    # get the testing target variable corresponding to the testing feature dates
    if get_end_of_week_for_rolling_average(opt.d_var) == 'Tue':
        # added on Mar 18: data_x_test.index are a single Friday observation
        data_y_test = data.set_index('date').loc[data_x_test.index, f'{opt.d_var}_t{opt.wk}_fwd']
    elif get_end_of_week_for_rolling_average(opt.d_var) == 'Thu':
        data_y_test = data.set_index('true_date').loc[data_x_test.index, f'{opt.d_var}_t{opt.wk}_fwd']
    
    # concatenate training features and target, drop any rows with missing values
    data_x_train, data_y_train = drop_na_xy(data_x_train, data_y_train)

    # return the prepared training and testing sets
    return data_x_train, data_y_train, data_x_test, data_y_test


def get_train_test_split_for_nontext(opt, X_train_start_date, X_train_end_date, x_test_date):
    # read in full dataset according to the specified dependent variable (d_var)
    data = data_set(opt.d_var)
    # create forward-looking dependent variables by shifting the specified weeks (wk) forward
    data[f'{opt.d_var}_t{opt.wk}_fwd'] = data[f'{opt.d_var}_t{opt.wk}'].shift(-opt.wk)
    if get_end_of_week_for_rolling_average(opt.d_var) == 'Tue':
        data_x = data[['date_Tue', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot', 'FutRet', 'StkIdx', 'xomRet', 'bpRet', 'rdsaRet',
          'tnote_10y', 'DFX', 'sp500Ret', 'basis', f'WIPI_{opt.wk}wk', 'trend', 'VIX', 'vix_diff', 'ovx_diff', 'RPsdf_growing', 'RPsdf_rolling',
          'BEME', 'Mom', 'BasMom', 'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']]\
                .set_index('date_Tue')
        # extract training data based on the start and end dates, sorting by date
        data_x_train = data_x[(data_x.index >= X_train_start_date) & (data_x.index <= X_train_end_date)]\
                        .sort_values('date_Tue')
    elif get_end_of_week_for_rolling_average(opt.d_var) == 'Thu':
        data_x = data[['true_date', 'DOilVol', 'OilVol', 'DInv', 'DProd', 'DSpot', 'FutRet', 'StkIdx', 'xomRet', 'bpRet', 'rdsaRet',
                  'tnote_10y', 'DFX', 'sp500Ret', 'basis', f'WIPI_{opt.wk}wk', 'trend', 'VIX', 'vix_diff', 'ovx_diff', 'RPsdf_growing', 'RPsdf_rolling',
                  'BEME', 'Mom', 'BasMom', 'DolBeta', 'InflaBeta', 'HedgPres', 'liquidity', 'OpenInt']]\
                        .set_index('true_date')
        # extract training data based on the start and end dates, sorting by date
        data_x_train = data_x[(data_x.index >= X_train_start_date) & (data_x.index <= X_train_end_date)]\
                        .sort_values('true_date')
    
    # added on Mar 18
    if get_end_of_week_for_rolling_average(opt.d_var) == 'Tue':
        data_x_test = data_x.loc[[x_test_date - pd.Timedelta('3 days')]]
        data_x_test.index += pd.Timedelta('3 days')
    elif get_end_of_week_for_rolling_average(opt.d_var) == 'Thu':
        data_x_test = data_x.loc[[x_test_date]]
    
    # added on Mar 11: if test row contains na, use the most recent observation
    data_x_train, data_x_test = fill_na_train_test(data_x_train, data_x_test)
    return data_x_train, data_x_test


def get_avg_univariate_r2(x_train, y_train, opt):
    kf = KFold(n_splits=opt.cvs)
    avg_r2_score = {}
    for d_var_one_predictor in x_train.columns:
        x_train_one_predictor = x_train[d_var_one_predictor]
        train_xy = pd.concat([x_train_one_predictor, y_train], axis=1).dropna()
        if len(train_xy) > len(x_train) * 0.9:
            x_train_one_predictor = train_xy.iloc[:,:-1]
            y_train = train_xy.iloc[:,-1]
            r2_scores = []
            for train_index, test_index in kf.split(x_train_one_predictor):
                X_cv_train, X_cv_test = x_train_one_predictor.iloc[train_index], x_train_one_predictor.iloc[test_index]
                X_cv_train = sm.add_constant(X_cv_train)
                X_cv_test = sm.add_constant(X_cv_test)
                y_cv_train, y_cv_test = y_train.iloc[train_index], y_train.iloc[test_index]
                model = lm.OLS(y_cv_train, X_cv_train, missing='drop').fit()
                y_cv_pred = model.predict(X_cv_test)
                r2_scores.append(r2_score(y_cv_test, y_cv_pred))
            avg_r2_score[d_var_one_predictor] = np.mean(r2_scores)
        else:
            avg_r2_score[d_var_one_predictor] = -999
    return avg_r2_score

def train_forward_lasso(x_train_text, x_train_nontext, y_train, opt):
    if opt.select_significant:
        avg_r2_score_text = get_avg_univariate_r2(x_train_text, y_train, opt)
        significant_text_vars = [key for key, value in sorted(avg_r2_score_text.items(), key=lambda item: item[1], reverse=True)[:opt.num_text_var]]
        
        avg_r2_score_nontext = get_avg_univariate_r2(x_train_nontext, y_train, opt)
        significant_nontext_vars = [key for key, value in sorted(avg_r2_score_nontext.items(), key=lambda item: item[1], reverse=True)[:opt.num_nontext_var]]
        
        significant_ind_vars = significant_text_vars + significant_nontext_vars
        x_train = pd.concat([x_train_text, x_train_nontext], axis=1)[significant_ind_vars]
        
    else:
        x_train = pd.concat([x_train_text, x_train_nontext], axis=1)
        
    x_train, y_train = drop_na_xy(x_train, y_train)

    ## Set up Lasso instance and grid search for penalty coefficient
    pre_model = Lasso(random_state=seed, fit_intercept=not opt.normalize)
    if opt.normalize:
        param_grid = [{'alpha': np.exp(np.linspace(-7, 0, 40))}]
    else:
        param_grid = [{'alpha': np.linspace(1e-5, 2, 40)}]
    grid_search = GridSearchCV(pre_model, param_grid, cv=opt.cvs, scoring='neg_mean_squared_error', n_jobs=-1)
    grid_search.fit(x_train, y_train)
    
    return grid_search


if __name__ == '__main__':
    
    opt = parse_option()
    print(opt)
    
    directory_name = f"wk_{opt.wk}_window_{opt.window}_frequency_{opt.frequency}_cvs_{opt.cvs}_rolling_{opt.rolling}_select_significant_{opt.select_significant}_num_nontext_var_{opt.num_nontext_var}_num_text_var_{opt.num_text_var}_normalize_{opt.normalize}"
    directory_path = os.path.join(f'res_{opt.model_name}', directory_name)
    os.makedirs(directory_path, exist_ok=True)
    
    time_col = data_set(opt.d_var)['date'] 
    time_lower = time_col[0]
    # a list of dates
    test_week_list = [time for time in time_col if time>=time_lower+pd.Timedelta(str(7*opt.window*52)+'days')][::opt.frequency]
    assert test_week_list[0].weekday() == 4, "The dates are not Fridays"
    
    res = {}

    for ii, end_of_train in enumerate(tqdm(test_week_list)): #[876:]
        try:
            if get_end_of_week_for_rolling_average(opt.d_var) == 'Tue':
                end_of_train -= pd.Timedelta('2 days')
            x_train_text, y_train, x_test_text, y_test = get_train_test_split_for_text(opt, end_of_train)
            X_train_start_date = x_train_text.index.min()
            X_train_end_date = x_train_text.index.max()
            x_test_date = y_test.index[0]
            x_train_nontext, x_test_nontext = get_train_test_split_for_nontext(opt, X_train_start_date, X_train_end_date, x_test_date)
            
            x_train_text_normed = standardize(x_train_text)
            x_test_text_normed = (x_test_text - x_train_text.mean())/x_train_text.std()
            x_train_nontext_normed = standardize(x_train_nontext)
            x_test_nontext_normed = (x_test_nontext - x_train_nontext.mean())/x_train_nontext.std()
            y_train_normed = standardize(y_train)
            
            model = train_forward_lasso(x_train_text_normed, x_train_nontext_normed, y_train_normed, opt)
            significant_ind_vars = list(model.feature_names_in_)
            x_train_normed_selected = pd.concat([x_train_text_normed, x_train_nontext_normed], axis=1)[significant_ind_vars]
            x_test_normed_selected = pd.concat([x_test_text_normed, x_test_nontext_normed], axis=1)[significant_ind_vars]
            
            x_train_normed_selected, y_train_normed = drop_na_xy(x_train_normed_selected, y_train_normed)
            pred_train_normed = model.predict(x_train_normed_selected)
            pred_test = model.predict(x_test_normed_selected) * y_train.std() + y_train.mean()
            assert len(pred_test) == len(y_test) == 1, 'Prediction or truth has wrong size'
            
            res[str(end_of_train)[:10]] = {'d_var':opt.d_var,
                                           'best_lambda':model.best_params_['alpha'], 
                                           'significant_ind_vars':significant_ind_vars,
                                            'pred':pred_test[0], 'true':y_test.astype(float).values[0],
                                            'mean':y_train.mean(), 'std':y_train.std(),
                                            'coefs':list(model.best_estimator_.coef_),
                                            'in_sample_R2':r2_score(y_train_normed.values, pred_train_normed)}
            # generate buckets of IS R2, and compare the OOS R2 for each bucket
            # still need ex-ante rules
        except IndexError:
            break
    torch.save(res, f'{directory_path}/{opt.d_var}.pt')