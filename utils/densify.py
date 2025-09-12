import random
import numpy as np


global_scale=1

def add_new_edges(H,u,hmap,avail,kval,local_kval):

    candidates = sorted(hmap.items(), key=lambda item: item[1], reverse=True)
    counts = min(len(candidates), kval - local_kval)

    select = 0
    i = 0
    while select < counts and i < len(candidates):
        idx = candidates[i][0]
        if avail[idx] == 1:
            H.add_edge(u, idx, weight=hmap[idx])
            select += 1

        i += 1

    return H


def Densify_v0(G, H, top_nodes, kval):

    n = G.number_of_nodes()
    avail = np.zeros(n)
    for u in H.nodes():
        avail[u] = 1

    for u in top_nodes:
        local_kval = len(list(H.neighbors(u)))
        hmap = {}
        for _ in range(20):
            stop = 0
            idx = u
            sc = 1
            steps = 0
            while stop != 1 and steps < 10:
                ngbrs = list(G.neighbors(idx))
                new_idx = random.choice(ngbrs)

                if new_idx in hmap:
                    hmap[new_idx] += sc
                else:
                    hmap[idx] = sc

                idx = new_idx
                sc = sc * global_scale
                steps += 1

        H=add_new_edges(H, u, hmap, avail, kval, local_kval)


    return H