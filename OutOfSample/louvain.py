import numpy as np
import random

def get_clusters(seed, graph):
    random.seed(seed)
    return graph.community_multilevel(weights='cosine')

def get_topic_vec(cluster, word_array, word2weight, weight=True):
    assert word_array[0] == 'fuel' # original order of words
    num_words = len(word_array)
    mask = np.zeros(num_words, dtype=bool)
    mask[cluster] = True
    if not weight:
        return mask
    else:
        topic_vec = np.zeros(num_words)
        topic_vec[mask] = [word2weight[word] for word in word_array[cluster]]
        return topic_vec
    
def get_matching_score(matrix, num_pairs='min', return_idx=False):
    
    """
    return the mean of greedily-searched best matching pairs between old and new topic vectors 
    based on max cos-similarity
    """
    
    if num_pairs == 'min':
        num_pairs = min(matrix.shape)
    elif num_pairs == 'max':
        num_pairs = max(matrix.shape)
        
    max_values, max_value_indices = [], []
    for _ in range(num_pairs):
        max_value_idx = np.unravel_index(np.argmax(matrix), matrix.shape)
        max_value = matrix[max_value_idx]
        max_values.append(max_value)
        matrix[max_value_idx[0], :] = -np.inf
        matrix[:, max_value_idx[1]] = -np.inf
        if return_idx:
            max_value_indices.append(max_value_idx)

    return max_values, max_value_indices
