import hnswlib
import numba
import numpy as np

def get_kNN(X, q=15):
    """
    Generate a k-nearest neighbors graph from the input data.
    :param X: Input data (numpy array).
    :param q: Number of nearest neighbors.
    :return: k-nearest neighbors list and distances.
    """
    n = X.shape[0]
    dim = X.shape[1]
    p = hnswlib.Index(space='l2', dim=dim)
    p.init_index(max_elements=n, ef_construction=200, M=64)
    p.add_items(X)
    p.set_ef(2*q)

    labels, dists = p.knn_query(X, k=q+1)
    knn_list = labels[:, 1:]
    knn_dists = dists[:, 1:]

    return knn_list, knn_dists


import networkx as nx
def get_KNN_nx(X,q=15):

    knn_list,_=get_kNN(X,q)
    G = nx.DiGraph()
    for i in range(len(knn_list)):
        for j in knn_list[i]:
            G.add_edge(i, j, weight=1)


    return G