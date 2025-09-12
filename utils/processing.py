import networkx as nx


def typename(o):
    return o.__class__.__name__

def children_of(o):
    return getattr(o, "children", []) or []


def cores(corespect):
    return [c for c in children_of(corespect) if typename(c) == "Core"]

def clustered_children(core):
    return [c for c in children_of(core) if typename(c) == "Clustered_Core"]

def propagated_children(clustered_core):
    return [p for p in children_of(clustered_core) if typename(p) == "Propagated_Data"]


def corespect_tracer(corespect):

    cs_tree=nx.DiGraph()

    for co in cores(corespect):
        cs_tree.add_edge(corespect, co)
        for cc in clustered_children(co):
            cs_tree.add_edge(co, cc)
            for pd in propagated_children(cc):
                cs_tree.add_edge(cc, pd)

    return cs_tree


def get_leaves(corespect):
    cs_tree=corespect_tracer(corespect)

    leaves = [node for node in cs_tree.nodes() if cs_tree.out_degree(node) == 0]

    return leaves


def _find_adaptive_start_point(arr):

    arr_set = set(arr)
    m = max(arr)

    skip_list=[x for x in range(min(arr), m + 1) if x not in arr_set]

    if len(skip_list) == 0:
        return -1

    elif len(skip_list) == 1:
        return skip_list[0]

    else:
        raise ValueError("Multiple skips found in round_info, which is unexpected.")


def extract_layers_labels(round_info,final_labels,mode='three_steps'):

    n=len(round_info)
    skip_pos = _find_adaptive_start_point(round_info)

    if skip_pos == -1:
        skip_pos = max(round_info) + 1

    if mode is None:
        raise ValueError("Mode must be specified.")

    elif mode=='three_steps':

        layer_num = 3
        layers = [[] for _ in range(layer_num)]
        labels = [[] for _ in range(layer_num)]

        for ell in range(n):

            if round_info[ell] == 0:
                layers[0].append(ell)
                labels[0].append(final_labels[ell])

            elif 0< round_info[ell] < skip_pos:
                layers[1].append(ell)
                labels[1].append(final_labels[ell])

            elif round_info[ell]>1:
                layers[2].append(ell)
                labels[2].append(final_labels[ell])

            else:
                continue

    elif mode=='layerwise':

        layer_num = len(set(round_info))-1 if -1 in round_info else len(set(round_info))
        layers = [[] for _ in range(layer_num)]
        labels = [[] for _ in range(layer_num)]

        for ell in range(n):

            if round_info[ell] == 0:
                layers[0].append(ell)
                labels[0].append(final_labels[ell])

            elif round_info[ell] ==-1:
                continue
                # layers[layer_num-1].append(ell)
                # labels[layer_num-1].append(final_labels[[ell]])

            elif round_info[ell] < skip_pos:
                layers[round_info[ell]].append(ell)
                labels[round_info[ell]].append(final_labels[ell])

            else:
                layers[round_info[ell]-1].append(ell)
                labels[round_info[ell]-1].append(final_labels[ell])


    else:
        raise ValueError("Mode not recognized.")

    return layers,labels


