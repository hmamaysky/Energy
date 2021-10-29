#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 20:16:51 2020

@author: billwu

This file plots the wordclouds of all the LDA topics in the paper in one graph.
"""
# %% Import Packages
import pandas as pd
import os
import wordcloud
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
plt.rcParams["font.family"] = "Times New Roman"

# %% Defining Functions
def topic_freq_dict(words):
    '''
    Take in the word-freq-topic dict, separate into topics and into word-freq dict
    '''
    # keys for the result dict
    topic_name_dict={'1':'Co','2':'Gom','3':'Env','4':'Epg',
                     '5':'Bbl','6':'Rpc','7':'Ep'}
    topic_topfreqs={}
    # for each topic, get word-freq dict and save as the value of the output dict
    for i in range(7):
        df_temp = words[words['Topic']==i+1].loc[:,['word','freq']].sort_values(by='freq',ascending=False)
        topic_topfreqs[topic_name_dict[str(i+1)]]=dict(zip(df_temp['word'],df_temp['freq']))
    return topic_topfreqs

def topic_wordCloud(freq_dict):
    '''
    Plot wordclouds of 7 topics in one graph
    '''
    # Arrange the plotting order of topics
    topic_name_dict={'1':'Co','2':'Gom','3':'Env','4':'Epg',
                     '5':'Bbl','6':'Rpc','7':'Ep'}
    # Get title of each topic
    topic_title_dict={'Co':'Company (Co)', 'Gom':'Global oil market (Gom)', 'Env':'Environment (Env)',
                      'Bbl':'Crude oil physical (Bbl)','Rpc':'Refining & petrochemicals (Rpc)',
                      'Ep':'Exploration & production (Ep)','Epg':'Energy/power generation (Epg)'}
    fig = plt.figure(figsize=(24,12),dpi=100)
    gs = gridspec.GridSpec(2,8,figure=fig)
    gs.update(wspace=0.3, hspace=0.5)
    plot_positions={1:[0,'0,2'], 2:[0,'2,4'], 3:[0,'4,6'], 4:[0,'6,8'],
                    5:[1,'1,3'], 6:[1,'3,5'], 7:[1,'5,7']}
    subplots={}
#    fig,axes=plt.subplots(4,2, figsize=(12,20),dpi=200)
    # i,j to adjust plotting position
#    i=0;j=0
    for k in range(7):
#        i=int(k/2)
#        j=k%2
        # get topic name
        topic = topic_name_dict[str(k+1)]
        # get word-freq dictionary for topic
        word_freq_lists=freq_dict[topic]
        # get topic title
        topic_title=topic_title_dict[topic]
        # instantiate the wc object
        wc=wordcloud.WordCloud(height=600, width=600, max_font_size=250, background_color='white')
        # fit the dict of current topic (word-freq)
        wc.fit_words(word_freq_lists)
        sub_p = plot_positions[k+1]
        subplots[k] = plt.subplot(gs[sub_p[0], int(sub_p[1].split(',')[0]):int(sub_p[1].split(',')[1])])
        subplots[k].imshow(wc, interpolation='bilinear')
        # turn off the axis
        subplots[k].axis('off')
        # set title and take care of fontsize, color as well as padding
        subplots[k].set_title(topic_title,fontsize=20,color='b',pad=35)
    # do not forget disguise the last plot
#    axes[3,1].axis('off')

    # adjust the position between subplots
#    plt.subplots_adjust(wspace=0.7,hspace=0.4)
    # save fig with some paddings at edges
    plt.savefig('../../../Paper/CEBRA 2020 PPT Figures/figures/Figure1.png', bbox_inches = 'tight',
                pad_inches = 0.5)
    
# %% Main Function
def main():
    wkdir = '/Users/billwu/Desktop/2020 Spring/2020 RA/Prof. Harry Mamaysky/'
    wkdir += 'Energy Project/Data Processing/LDA analysis/Louvain method'
    os.chdir(wkdir)
    # Read in the df for word-freq-topic
    words = pd.read_csv('clustering_C_202007.csv')
    # Separate df into 7 word-freq dicts according to their own topic
    topic_freqdicts = topic_freq_dict(words)
    # Plot and save the figure
    topic_wordCloud(topic_freqdicts)
    
# %% Main Function
if __name__=='__main__':
    main()