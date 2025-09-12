
from ..utils.processing import extract_layers_labels
from ..utils.processing import typename
from ..utils import processing as tracer

from .deg_analysis import deg_one_layer
from .deg_analysis import  de_two_groups
from .deg_analysis import calc_deg
import numpy as np

from ..utils.processing import get_leaves


class Deg:
    def __init__(
            self,
            corespect=None,
            X_cg=None,
            gene_idx=None,
            true_labels=None,
            all_gene_names=None,
            top_deg_num=10,
            mode='three_steps'
    ):

        self.corespect = corespect
        self.X_cg = X_cg
        self.gene_idx = gene_idx
        self.true_labels = true_labels
        self.all_gene_names = all_gene_names
        self.top_deg_num = top_deg_num
        self.mode = mode




    def deg_two_sets(self,prop_data, cluster1,cluster2,top_deg_num_per_layer=None,custom_labels=None):


        if top_deg_num_per_layer is None:
            top_deg_num_per_layer=self.top_deg_num

        if typename(prop_data) != "Propagated_Data":
            raise ValueError("Please provide a Propagated_Data object.")

        if self.X_cg is None:
            raise ValueError("Please provide the cell-gene matrix X_cg.")

        if self.gene_idx is None:
            self.gene_idx = [i for i in range(self.X_cg.shape[1])]


        if custom_labels is None:
            final_labels=prop_data.final_labels
        else:
            final_labels=custom_labels

        round_info=prop_data.round_info


        if set(cluster1).intersection(set(final_labels)) ==0 or set(cluster2).intersection(set(final_labels))==0:

            raise ValueError("Cluster1 not found in the final labels.")



        layers, labels = extract_layers_labels(round_info, final_labels, mode=self.mode)


        points_full = []
        gene_idx=[]


        #Get the complete list of all genes that are differentially expressed in any layer (top deg_num genes in each layer)
        for layer_num, points_round in enumerate(layers, start=0):

            points_full.extend(points_round)

            labels_set=set(labels[layer_num])

            if set(cluster1) & labels_set and set(cluster2) & labels_set:



                data1=de_two_groups(points_full,final_labels,self.gene_idx, self.X_cg,cluster1,cluster2,top_deg_num_per_layer=top_deg_num_per_layer)
                data2=de_two_groups(points_round,final_labels,self.gene_idx,self.X_cg,cluster1,cluster2,top_deg_num_per_layer=top_deg_num_per_layer)

                for gidx in data1.keys():
                    gene_idx.append(gidx)

                for gidx in data2.keys():
                    gene_idx.append(gidx)



        print("Total candidate genes found (with possible duplicates): ",len(gene_idx))
        gidx_set=np.array(list(set(gene_idx))).astype(int)
        print(len(gidx_set)," unique genes found in total.")
        print("\n")


        #Now we obtain the final dictionary.
        final_data={}
        for idx in gidx_set:
            final_data[idx]= {}


        points_full = []
        for layer_num, points_round in enumerate(layers, start=0):

            points_full.extend(points_round)


            labels_set=set(labels[layer_num])

            if set(cluster1) & labels_set and set(cluster2) & labels_set:

                temp_data=de_two_groups(points_full, final_labels, gidx_set, self.X_cg, cluster1, cluster2, top_deg_num_per_layer=len(gidx_set))

                for idx in temp_data.keys():
                    final_data[idx][str(layer_num)+'+']=temp_data[idx]

                temp_data = de_two_groups(points_round, final_labels, gidx_set, self.X_cg, cluster1, cluster2,
                                          top_deg_num_per_layer=len(gidx_set))

                for idx in temp_data.keys():
                    final_data[idx][str(layer_num)] = temp_data[idx]




        return final_data




    def detailed_deg(self,cluster1,cluster2,top_deg_num_per_layer=None,custom_labels=None):

        pd_list=get_leaves(self.corespect)

        deg_list=[]

        for prop_data in pd_list:

            if custom_labels is None:
                print("Processing a Propagated_Data object with final labels: ",set(prop_data.final_labels))
            else:
                print("Processing a Propagated_Data object with custom labels: ",set(custom_labels))

            temp_dict=self.deg_two_sets(prop_data, cluster1,cluster2,top_deg_num_per_layer=top_deg_num_per_layer,custom_labels=custom_labels)
            deg_list.append(temp_dict)


        return deg_list



