

from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

import operator
import numpy as np

def calc_deg(set1, set2, gene_idx_list, matrix):

    stat_dict=[]

    print(len(set1),len(set2),type(set1),type(set2))

    x=0
    for g in gene_idx_list:
        x=x+1

        group1=matrix[set1,g]
        group2=matrix[set2,g]




        mean1=group1.mean()
        mean2=group2.mean()

        pct1=np.count_nonzero(group1)/len(group1)
        pct2=np.count_nonzero(group2)/len(group2)

        if mean1>0.1 and mean1-mean2>0.05:

            if mean2==0:
                mean2=1/len(set2)


            stat, p_value = ranksums(group1, group2)
            logfc=np.log2(mean1/mean2)
            stat_dict.append([g,stat,p_value,logfc,group1.mean(),group2.mean(),pct1,pct2])


        if x%10000==0:
            print(x)


    if len(stat_dict)==0 :
       print("No suitable DEG found")
       return None


    stat_dict=np.array(stat_dict)
    p_values=stat_dict[:,2]


    rej,p_adj,_,_ = multipletests(p_values, method='fdr_bh')

    # print(rej)
    # print(len(stat_dict[:,2]),len(p_adj))

    stat_dict[:,2]=p_adj

    return stat_dict


#points are a subset of the n points
#label is a dictionary of (point,label of point) format

import matplotlib.pyplot as plt
from tabulate import tabulate


def deg_one_layer(points, label, gene_idx, matrix, top_deg_num, all_gene_names=None):

    layer_deg_dict={}

    if all_gene_names is None:
        all_gene_names=["G"+str(i) for i in range(matrix.shape[1])]

    c_num=len(set(label))

    label_set=[[] for _ in range(c_num)]

    label_hmap={}
    rev_label_hmap={}
    t=0
    for ell in set(label):
        label_hmap[ell]=t
        rev_label_hmap[t]=ell
        t=t+1

    label_num=[]
    for ell in range(len(label)):
        label_num.append(label_hmap[label[ell]])


    for i in points:
        idx=label[i]
        label_set[label_hmap[idx]].append(i)



    for ell in set(label)-{-1} if -1 in set(label) else set(label):

        print("Class in consideration=",ell)

        idx=label_hmap[ell]

        set1=label_set[idx]
        set2=list(set(points) - set(set1))

        print("size of class and rest",len(set1),len(set2))



        stat_dict=calc_deg(set1, set2, gene_idx, matrix)

        if stat_dict is None:
            print("No DEG found, skip")
            continue

        #changed from 1 to 3 here.
        sorted_dict = sorted(stat_dict, key=operator.itemgetter(3), reverse=True)


        # fig, axes = plt.subplots(int(np.ceil(top_deg_num / 5)), 5, figsize=(20, 8))
        # # Flatten the axes array for easy iteration
        # axes = axes.flatten()



        headers=["gene#"," gene"," zs"," pval","logfc","mean1","mean2"]

        data=[headers]


        for i in range(min(top_deg_num,len(sorted_dict))):

            gidx,zs,pval,logfc,mean1,mean2=list(sorted_dict[i])
            #print(f" {gidx:.1f}, {all_gene_names[int(gidx)]},{zs:.2f},{pval:.3f},{logfc:.3f},{mean1:.3f},{mean2:.3f}")

            data.append([gidx,all_gene_names[int(gidx)],zs,pval,logfc,mean1,mean2])

            z=matrix[set1,int(gidx)]
            lz=np.array(label_num)[set1]
            y=matrix[set2,int(gidx)]
            ly=np.array(label_num)[set2]

            xx=[i for i in range(len(points))]
            yx=list(z)+list(y)
            lf=list(lz)+list(ly)

            # axes[i].scatter(xx,yx,c=lf,s=2)
            # axes[i].set_title(str(ell)+' '+all_gene_names[int(gidx)])

        layer_deg_dict[ell]=data

#        print(tabulate(data, headers="firstrow", tablefmt="grid",floatfmt=".3f"))

        #plt.tight_layout()
        #plt.show()

#        print("\n\n\n\n\n")

    return layer_deg_dict



#cluster 1 is a list and cluster 2 is a list
def de_two_groups(points,final_labels,gene_idx, matrix,cluster1,cluster2,top_deg_num_per_layer=None):

    if top_deg_num_per_layer is None:
        top_deg_num_per_layer=10



    set1=[]
    set2=[]

    for ell in points:
        if final_labels[ell] in cluster1:
            set1.append(ell)
        elif final_labels[ell] in cluster2:
            set2.append(ell)


    stat_dict = calc_deg(set1, set2, gene_idx, matrix)


    if stat_dict is None:
        print("No DEG found, skip")
        return {}


    sorted_dict = sorted(stat_dict, key=operator.itemgetter(3), reverse=True)




    data = {}
    for i in range(min(len(sorted_dict), top_deg_num_per_layer)):
        gidx, zs, pval, logfc, mean1, mean2,pct1,pct2 = list(sorted_dict[i])
        data[gidx]=[pval,logfc,mean1,mean2,pct1,pct2]



    return data



