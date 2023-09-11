#!/user/kh3191/.conda/envs/nlp/bin/python

import igraph as ig
import pandas as pd
import random
import numpy as np
import torch
import os
from louvain import get_clusters
from cosine import get_cosine, get_word_set
from deap import base, creator, tools, algorithms


from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

seed = 12345
random.seed(seed)
np.random.seed(seed)

import argparse
from argparse import RawTextHelpFormatter
def parse_option():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--inputWordsPath', type=str, 
           default='clustering_C.csv')
    parser.add_argument('--dtmPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/dtm_numeric_441')
    parser.add_argument('--outputPath', type=str, 
           default='/shared/share_mamaysky-glasserman/energy_drivers/2023/DataProcessing/article_measure/rolling_topic')
    parser.add_argument('--updating_window_in_years', type=int, default=5)
    parser.add_argument('--min_topic_words', type=int, default=10)
    parser.add_argument('--remainder_topic', type=int, default=999)
    parser.add_argument('--GA_npop', type=int, default=1000)
    parser.add_argument('--GA_ngen', type=int, default=50)
    parser.add_argument('--GA_prob_crossing', type=float, default=0.5)
    parser.add_argument('--GA_prob_mutating', type=float, default=0.2)
    parser.add_argument('--GA_record_stats', type=bool, default=False)
    parser.add_argument('--rolling_index', type=int, default=0)
    opt = parser.parse_args()
    return opt
    
def allocate_small_topics_genetic(graph, membership_new, remainder_topic,
                                  npop, ngen, prob_crossing, prob_mutating,
                                  record_stats):
    
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    toolbox = base.Toolbox()

    existing_topic = list(set(membership_new).difference([remainder_topic]))
    def init_individual():
        return [random.choice(existing_topic) if x == remainder_topic else x for x in membership_new]

    toolbox.register("individual", tools.initIterate, creator.Individual, init_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    def evaluate(individual):
        return graph.modularity(individual, weights='cosine'),
    
    def mutate(individual, mutate_prob=0.1):
        for i in range(len(individual)):
            if random.random() < mutate_prob:
                individual[i] = random.choice(existing_topic)
        return individual,

    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutate)
    toolbox.register("select", tools.selTournament, tournsize=3)
    toolbox.register("evaluate", evaluate)

    if record_stats:
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("min", min)
        stats.register("max", max)
        verbose = True
    else:
        stats = None
        verbose = None
    
    population = toolbox.population(npop)
    hall_of_fame = tools.HallOfFame(maxsize=1)
    algorithms.eaSimple(population, toolbox, prob_crossing, prob_mutating, ngen, stats=stats, 
                        halloffame=hall_of_fame, verbose=verbose)

    population.sort(key=lambda x: x.fitness.values[0], reverse=True)
    index_99th_percentile = int(len(population) * 0.01)
    best_individual = population[index_99th_percentile]
    best_modularity = evaluate(best_individual)[0]
    return best_individual, best_modularity


if __name__ == "__main__":
    opt = parse_option()
    print(opt)
    
    word_set = get_word_set(opt.inputWordsPath, silent=True)
    pathin = [f"{opt.dtmPath}/{f}" for f in os.listdir(opt.dtmPath)]
    pathin.sort()

    res = {}
    res['mod'], res['membership'] = {}, {}

    for update_window in [opt.updating_window_in_years]:
        max_rolling_index = len(pathin)-update_window*12
        #for i in tqdm(range(len(pathin)-update_window*12+1)):
        for i in [opt.rolling_index]:
            assert max_rolling_index == 266
            assert 0 <= i <= max_rolling_index
            YYYYMM = pathin[i+update_window*12-1][-14:-8]
#             for key in res.keys():
#                 res[key][YYYYMM] = []
                
            frames = [pd.read_csv(j, delimiter=',') for j in pathin[i:i+update_window*12]]
            df_rolling = pd.concat(frames)
            df_rolling = df_rolling.query('words>=0')
            df_cosine = get_cosine(df_rolling, word_set, silent=True)

            graph = ig.Graph.Weighted_Adjacency(df_cosine.values, mode='undirected', attr='cosine', loops=False)

            mod_list, membership_list = [], []
            for s in tqdm(range(10)):
                c = get_clusters(s, graph)
                filtered_idx = set(filter(lambda x: c.sizes()[x] < opt.min_topic_words, range(len(c.sizes()))))
                membership_new = list(map(lambda x: opt.remainder_topic if x in filtered_idx else x, c.membership))
                best_membership, best_modularity = allocate_small_topics_genetic(graph, 
                                                                                 membership_new, 
                                                                                 opt.remainder_topic, 
                                                                                 npop=opt.GA_npop, 
                                                                                 ngen=opt.GA_ngen, 
                                                                                 prob_crossing=opt.GA_prob_crossing,
                                                                                 prob_mutating=opt.GA_prob_mutating,
                                                                                 record_stats=opt.GA_record_stats)

                mod_list.append(best_modularity) 
                membership_list.append(best_membership)
            best_s = np.argmax(mod_list)
            res['mod'][YYYYMM] = mod_list[best_s]
            res['membership'][YYYYMM] = membership_list[best_s]

    torch.save(res, f'{opt.outputPath}/res_{YYYYMM}.pt')
