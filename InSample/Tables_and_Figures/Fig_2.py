#!/user/kh3191/.conda/envs/nlp/bin/python
"""
adapted from billwu

This file plots the freq, sent of textual vars, and some other related quantities like PCA, unusualness, and artcount, etc.
"""
# %% Import Packages
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import seaborn as sns
sns.set(rc={"axes.labelsize": 15, 
            "axes.labelweight": 'bold',
            "axes.titlesize": 20,
            "axes.titleweight": 'bold',
            "xtick.labelsize": 15, 
            "ytick.labelsize": 15})
# plt.rcParams["font.family"] = "Times New Roman"
from datetime import datetime, timedelta

def generate_date_list(start, end):
    if isinstance(start, int):
        start_date = datetime.strptime(str(start), '%Y%m%d')
    else:
        start_date = datetime.strptime(start, '%Y-%m-%d')
    if isinstance(end, int):
        end_date = datetime.strptime(str(end), '%Y%m%d')
    else:
        end_date = datetime.strptime(end, '%Y-%m-%d')
    date_list = []
    while start_date <= end_date:
        date_list.append(int(start_date.strftime('%Y%m%d')))
        start_date += timedelta(days=1)
    return date_list

# need to modify the text: topic > 0.25
def get_headlines(df_all, topic, event_dates, use_list=False, top=None):
    start_date, end_date = event_dates[str(topic)][0]
    if use_list:
        date_list = generate_date_list(start_date, end_date)
        filtered_df = df_all.query(
            f"`date` in {date_list} and `entropy` >= 2 and `total` >= 100 and `Topic{topic}` > 0.25"
        ).copy()
    else:
        if isinstance(start_date, str):
            start_date = int(datetime.strptime(start_date, '%Y-%m-%d').strftime('%Y%m%d'))
        if isinstance(end_date, str):
            end_date = int(datetime.strptime(end_date, '%Y-%m-%d').strftime('%Y%m%d'))
        filtered_df = df_all.query(
            f"@start_date <= date <= @end_date and `entropy` >= 2 and `total` >= 100 and `Topic{topic}` > 0.25"
        ).copy()
    filtered_df['weight'] = filtered_df[f'Topic{topic}'] * filtered_df['sentiment']
    filtered_df = filtered_df.sort_values('weight', ascending=True)
    filtered_df = filtered_df[['sentiment', 'entropy', 'date', 'headline']]
    if top:
        filtered_df = filtered_df.head(top)
    return filtered_df

def format_float(value):
    if isinstance(value, float):
        return '{:.3f}'.format(value)
    else:
        return value
    
def int_to_date(value):
    if isinstance(value, int):
        date_obj = datetime.strptime(str(value), '%Y%m%d')
        return date_obj.strftime('%-m/%-d/%Y')
    else:
        return value

def find_next_wednesday(date_str, date_format="%Y-%m-%d"):
    given_date = datetime.strptime(date_str, date_format)
    weekday = given_date.weekday()
    days_to_add = (2 - weekday) % 7  # 2 corresponds to Wednesday
    next_wednesday = given_date + timedelta(days=days_to_add)
    return next_wednesday.strftime(date_format)

def find_previous_wednesday(date_str, date_format="%Y-%m-%d"):
    given_date = datetime.strptime(date_str, date_format)
    weekday = given_date.weekday()
    days_to_subtract = (weekday - 2) % 7  # 2 corresponds to Wednesday
    previous_wednesday = given_date - timedelta(days=days_to_subtract)
    return previous_wednesday.strftime(date_format)

# Arrange the plotting order of topics
topic_name_dict={'1':'Co','2':'Gom','3':'Env','4':'Epg',
                 '5':'Bbl','6':'Rpc','7':'Ep'}
# Get title of each topic
# topic_title_dict={'Co':'Company (Co)', 'Gom':'Global Oil Market (Gom)', 'Env':'Environment (Env)',
#                   'Bbl':'Crude Oil Physical (Bbl)','Rpc':'Refining & Pertrochemicals (Rpc)',
#                   'Ep':'Exploration & Production (Ep)','Epg':'Energy/Power Generation (Epg)'}
# old topic names
old2new = {'Env':'Dist', 'Epg':'Ge'}
topic_title_dict={'Co':'Company (Co)', 
                  'Gom':'Global Oil Market (Gom)', 
                  'Dist':'Distribution (Dist)',
                  'Bbl':'Crude Oil Physical (Bbl)',
                  'Rpc':'Refining & Pertrochemicals (Rpc)',
                  'Ep':'Exploration & Production (Ep)',
                  'Ge':'Generation & Environment (Ge)'}

# %% Defining Functions
def plot_freq(dataset, event_dates, yscale=1, rectangle=False):
    # ytick dict by topic
#     topic_ytick_dict={'Co':[0,0.1,0.2,0.3,0.4], 'Gom':[0,0.1,0.2,0.3,0.4,0.5,0.6], 'Env':[0,0.05,0.1,0.15],
#                       'Bbl':[0,0.05,0.1],'Rpc':[0,0.01,0.02,0.03,0.04],
#                       'Ep':[0,0.05,0.1,0.15],'Epg':[0.2,0.3,0.4,0.5]}
    fig, axes = plt.subplots(4,2,figsize=(16,20),dpi=200)
    for k in range(7):
        i = int(k/2)
        j = k%2
        # get topic name
        topic = topic_name_dict[str(k+1)]
        # get topic title
        try:
            topic_title = topic_title_dict[topic]
        except KeyError:
            topic_title = topic_title_dict[old2new[topic]]
        
        axes[i,j].plot(dataset['date'], dataset['ftopic'+str(k+1)+'_4wk'], color='b')
        # turn off the axis
        # set title and take care of color as well as padding
        axes[i,j].set_title(topic_title, 
                            color='k', 
                            weight='bold', 
                            backgroundcolor='#f0f0f0', 
                            pad=15)
        axes[i,j].set_xlim(dataset['date'].values[0],dataset['date'].values[-1])
        
        if event_dates[str(k+1)]:
            event_date = event_dates[str(k+1)][0]
            if isinstance(event_date, str):
                axes[i,j].plot(datetime.strptime(event_date,'%Y-%m-%d'),
                    yscale*dataset.loc[dataset['date_Wed'] == event_date,'ftopic'+str(k+1)+'_4wk'],
                    marker='*',color='red',markersize=16)
            else:
                start_date, end_date = event_date
                highlight_dataset = dataset.query('@start_date <= date_Wed <= @end_date').copy()
                axes[i,j].plot(highlight_dataset['date_Wed'],
                        yscale*highlight_dataset['ftopic'+str(k+1)+'_4wk'],
                        color='red', linewidth=3)
                for date in [find_next_wednesday(start_date), find_previous_wednesday(end_date)]:
                    axes[i,j].plot(datetime.strptime(date,'%Y-%m-%d'),
                            yscale*dataset.loc[dataset['date_Wed'] == date,'ftopic'+str(k+1)+'_4wk'],
                            marker='*',color='red',markersize=16)

    # do not forget disguise the last plot
    axes[3,1].axis('off')
    
    if rectangle:
        rect = Rectangle((0, 0), 1, 1, linewidth=3, edgecolor='k', facecolor='none', transform=fig.transFigure, clip_on=False)
        fig.patches.extend([rect])
    fig.suptitle('Panel A: Topical Frequency', fontsize=30, y=1)
    fig.subplots_adjust(top=0.90)
    # adjust the position between subplots
    fig.subplots_adjust(wspace=0.2,hspace=0.4)
    fig.tight_layout()
    fig.savefig('Fig2_1.pdf')
    
def plot_sent(dataset, event_dates, yscale=1000, rectangle=False):
    # y position of the scale annotation (\times 10^{-3})
#     annotation_y_dict={'Co':0.1, 'Gom':0.24, 'Env':0.18,
#                       'Bbl':0.05,'Rpc':0.06,
#                       'Ep':0.1,'Epg':0.26}
    fig, axes = plt.subplots(4,2,figsize=(16,20),dpi=200)
    for k in range(7):
        i = int(k/2)
        j = k%2
        # get topic name
        topic = topic_name_dict[str(k+1)]
        # get topic title
        try:
            topic_title = topic_title_dict[topic]
        except KeyError:
            topic_title = topic_title_dict[old2new[topic]]
        
        axes[i,j].plot(dataset['date'], yscale*dataset['stopic'+str(k+1)+'_4wk'], color='b')
        # turn off the axis
        # set title and take care of color as well as padding
        axes[i,j].set_title(topic_title, 
                            color='k', 
                            weight='bold', 
                            backgroundcolor='#f0f0f0', 
                            pad=15)
        axes[i,j].set_xlim(dataset['date'].values[0],dataset['date'].values[-1])
        ymin, ymax = axes[i,j].get_ylim()
        axes[i,j].text(pd.Timestamp('1998-05-01'), ymin, r'$\times 10^{-3}$', fontsize=12)
        #axes[i,j].set_ylabel('x$10^{-3}$', rotation=0, labelpad=20, verticalalignment='center', fontsize=12)
        
        if event_dates[str(k+1)]:
            event_date = event_dates[str(k+1)][0]
            if isinstance(event_date, str):
                axes[i,j].plot(datetime.strptime(event_date,'%Y-%m-%d'),
                    yscale*dataset.loc[dataset['date_Wed'] == event_date,'stopic'+str(k+1)+'_4wk'],
                    marker='*',color='red',markersize=16)
            else:
                start_date, end_date = event_date
                highlight_dataset = dataset.query('@start_date <= date_Wed <= @end_date').copy()
                axes[i,j].plot(highlight_dataset['date_Wed'],
                        yscale*highlight_dataset['stopic'+str(k+1)+'_4wk'],
                        color='red', linewidth=3)
                for date in [find_next_wednesday(start_date), find_previous_wednesday(end_date)]:
                    axes[i,j].plot(datetime.strptime(date,'%Y-%m-%d'),
                            yscale*dataset.loc[dataset['date_Wed'] == date,'stopic'+str(k+1)+'_4wk'],
                            marker='*',color='red',markersize=16)
                    
    # do not forget disguise the last plot
    axes[3,1].axis('off')
    
    if rectangle:
        rect = Rectangle((0, 0), 1, 1, linewidth=3, edgecolor='k', facecolor='none', transform=fig.transFigure, clip_on=False)
        fig.patches.extend([rect])
    fig.suptitle('Panel B: Topical Sentiment', fontsize=30, y=1)
    fig.subplots_adjust(top=0.90)
    # adjust the position between subplots
    fig.subplots_adjust(wspace=0.2,hspace=0.4)
    fig.tight_layout()
    fig.savefig('Fig2_2.pdf')

def plot_others(dataset, rectangle=False):
    # Arrange the plotting order of topics
    topic_name_dict={'1':'artcount_4wk','2':'entropy_4wk','3':'PCAsent','4':'PCAfreq',
                     '5':'PCAall'}
    # Get title of each topic
    topic_title_dict={'artcount_4wk':'Article Counts', 'entropy_4wk':'Entropy', 'PCAsent':'First PC of Normalized Topical Sentimant',
                      'PCAfreq':'First PC of Normalized Topical Frequency','PCAall':'First PC of Normalized all Text Variables'}
    fig, axes = plt.subplots(4,2,figsize=(16,20),dpi=200)
    for k in range(5):
        i = int(k/2)
        j = k%2
        # get topic name
        topic = topic_name_dict[str(k+1)]
        # get topic title
        topic_title = topic_title_dict[topic]
        
        axes[i,j].plot(dataset['date'], dataset[topic_name_dict[str(k+1)]], color='b')
        # turn off the axis
        # set title and take care of color as well as padding
        axes[i,j].set_title(topic_title, 
                            color='k', 
                            weight='bold', 
                            backgroundcolor='#f0f0f0', 
                            pad=15)
        axes[i,j].set_xlim(dataset['date'].values[0],dataset['date'].values[-1])
        
    # do not forget disguise the last plot
    axes[2,1].axis('off')
    axes[3,0].axis('off')
    axes[3,1].axis('off')
    
    if rectangle:
        rect = Rectangle((0, 0), 1, 1, linewidth=3, edgecolor='k', facecolor='none', transform=fig.transFigure, clip_on=False)
        fig.patches.extend([rect])
    fig.suptitle('Panel C: Article Counts, Unusualness and PCA series', fontsize=30, y=1)
    fig.subplots_adjust(top=0.90)
    # adjust the position between subplots
    fig.subplots_adjust(wspace=0.2,hspace=0.4)
    fig.tight_layout()
    fig.savefig('Fig2_3.pdf')

def get_dataset(file_dir):
    # read the latest dataset
    dataset = pd.read_stata(file_dir)
    date_cols_price = [x for x in list(dataset.columns.values) if 'date' in x]
    dataset = dataset.rename(columns={x:'_'.join(x.split('_')[:-1]) for x in set(dataset.columns.values) if x not in date_cols_price})

    # remedy missing var in in-sample analysis
    dataset = dataset.rename(columns={'date_Fri':'date'})
    dataset['sent'] = dataset['sCo']+dataset['sGom']+dataset['sEnv']\
                     +dataset['sEpg']+dataset['sBbl']+dataset['sRpc']+dataset['sEp']
    return dataset

# %% Main Function
def main():

    event_dates = {'1':('2000-09-20','UK fuel protests'),
               '2':('2002-04-24','Failed Venezuelan coup'),
               '3':('2015-10-14','Volkswagen emissions scandal'),
               '4':('2002-02-13','Post-bankruptcy Enron hearings'),
               '5':('2005-09-21','Hurricane Katrina'),
               '6':None,
               '7':('2010-06-02','BP oil spill aftermath')}
    dataset = get_dataset('transformed_data_prices_v19.dta')
    plot_freq(dataset, event_dates)
    plot_sent(dataset, event_dates)
    plot_others(dataset)
    
# %% Main Process
if __name__=='__main__':
    main()