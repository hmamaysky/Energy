#!/user/kh3191/.conda/envs/nlp/bin/python
import pandas as pd
import numpy as np
seed = 12345
np.random.seed(seed)
import statsmodels.api as sm
import statsmodels.regression.linear_model as lm
from sklearn.linear_model import Lasso
import os

from get_best_lambda import get_train_test_split
import torch

import re
from collections import Counter
from tqdm import tqdm

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
    parser.add_argument('--top_R2', type=int, default=2)
    parser.add_argument('--normalize', type=bool, default=True)
    opt = parser.parse_args()
    return opt

df_cluster = pd.read_csv('TopicLabel_Cluster.csv').sort_values('X')
df_cluster.index = df_cluster['X'].apply(lambda x: x[:4]+x[5:7]+x[-2:])
df_cluster = df_cluster['cluster']
all_features = ['artcount', 'entropy'] + [f'meta_f{i+1}' for i in range(df_cluster.max())]\
                + [f'meta_s{i+1}' for i in range(df_cluster.max())]\
                + ['PCAsent', 'PCAall', 'PCAfreq']

def get_meta_label(lasso_features):

    lasso_features_cluster = []
    for lasso_feature in lasso_features:
        if re.match(r'f\d+$', lasso_feature):
            lasso_features_cluster.append('meta_f'+str(df_cluster[f"{YYYYMM_end}_{int(lasso_feature[1:])-1}"]))
        elif re.match(r's\d+$', lasso_feature):
            lasso_features_cluster.append('meta_s'+str(df_cluster[f"{YYYYMM_end}_{int(lasso_feature[1:])-1}"]))
        else:
            lasso_features_cluster.append(lasso_feature)
            
    return lasso_features_cluster

if __name__ == '__main__':
    
    opt = parse_option()
    print(opt)
    
    directory_name = f"wk:{opt.wk}_window:{opt.window}_frequency:{opt.frequency}_cvs:{opt.cvs}_rolling:{opt.rolling}_select_significant:{opt.select_significant}_top_R2:{opt.top_R2}_normalize:{opt.normalize}"
    directory_path = os.path.join('res_Forward_Lasso', directory_name, 'test')
    os.makedirs(directory_path, exist_ok=True)

    d_var_list = ['FutRet', 'xomRet', 'bpRet', 'rdsaRet', 'DSpot', 'DOilVol', 'DInv', 'DProd']
    res = {}

    for d_var in d_var_list:
        
        print(d_var)
        res[d_var] = {'MSE': {}, 'true': {}, 'pred': {}, 'pred_wls': {}, 'trend': {}}

        res_train = torch.load(f"res_Forward_Lasso/{directory_name}/train/{d_var}.pt")
        best_lambda_dic = res_train['best_lambda_dic']
        significant_ind_vars_dic = res_train['significant_ind_vars_dic']

        frequency_all_features_list = []
        forecast_start_list = []

        for forecast_start, (YYYYMM_end, best_lambda, scaler) in tqdm(best_lambda_dic.items()):
            x_train, data_xtest, y_train, data_ytest, YYYYMM_end, scaler = get_train_test_split(d_var, pd.to_datetime(forecast_start), opt)
            
            if opt.select_significant:
                significant_ind_vars = significant_ind_vars_dic[forecast_start]
                x_train = x_train[significant_ind_vars]
                data_xtest = data_xtest[significant_ind_vars]
            x_test = data_xtest
            y_test = data_ytest[f'{d_var}_t{opt.wk}'].values

            if best_lambda:
                ## Update the coefficients using the selected penalty 
                reg = Lasso(alpha=best_lambda, random_state=seed, fit_intercept=not opt.normalize)
                reg.fit(x_train, y_train)
                selected_features = reg.feature_names_in_[reg.coef_ != 0]
                selected_features_cluster = get_meta_label(selected_features)

                # main parts for predictions
                res[d_var]['true'][forecast_start] = y_test.squeeze().astype(float)
                res[d_var]['pred'][forecast_start] = reg.predict(x_test).squeeze() * scaler['ystd'] + scaler['ymean']
                res[d_var]['trend'][forecast_start] = scaler['ymean']

                diff = reg.predict(x_test) - y_test
                res[d_var]['MSE'][forecast_start] = (diff ** 2).sum()

                model = lm.WLS(y_train,x_train, missing='drop',weights=np.linspace(1,2,len(x_train))).fit()
                res[d_var]['pred_wls'][forecast_start] = model.predict(x_test).squeeze() * scaler['ystd'] + scaler['ymean']


            frequency_all_features = {var: 0 for var in all_features}
            frequency_all_features.update(Counter(selected_features_cluster))
            frequency_all_features_list.append(list(frequency_all_features.values()))
            forecast_start_list.append(forecast_start)


        if opt.rolling:
            df_all_features = pd.DataFrame(frequency_all_features_list, columns=all_features, index=forecast_start_list).T
            res[d_var]['df_all_features'] = df_all_features
        torch.save(res, f"{directory_path}/res.pt")
