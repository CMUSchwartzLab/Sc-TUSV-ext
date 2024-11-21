import sys
import os
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist
from sklearn.metrics import adjusted_rand_score

def compute_l1_distances(C):
    return cdist(C, C, metric='cityblock') # manhattan/L1 distance

def compute_k_means(distances, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(distances)
    return kmeans.labels_

def get_true_labels(Z_file):
    Z = pd.read_csv(Z_file,header=None, sep='\t')
    return np.array([np.where(Z.iloc[i]==1)[0][0] for i in range(Z.shape[0])])

def save_clusters_file_sctusv(pred_labels, sampleList, out_dir):
    rows = []
    
    for i in range(len(pred_labels)):
        rows.append({'Cells': sampleList[i].split(".")[0], 'cluster': str(pred_labels[i]+1)})
    df = pd.DataFrame(rows)
    df.to_csv(out_dir+"/pred_kmeans_clusters.tsv",sep='\t',index=False)

def save_clusters_file(pred_labels, out_file_path):
    rows = []
    for i in range(len(pred_labels)):
        rows.append({'Cells': 'sample' + str(i+1), 'cluster': str(pred_labels[i]+1)})
    df = pd.DataFrame(rows)
    df.to_csv(out_file_path+"/pred_kmeans_clusters.tsv",sep='\t',index=False)
    
#print("Command to run: python cluster_l1.py input_C_file_path true_Z_file_path output_file_path n_clusters")
#C_file = sys.argv[1]
#Z_file = sys.argv[2]
#out_file = sys.argv[3]
#n_clusters = int(sys.argv[4])

# getting predicted labels and save a c2cl file
#C = pd.read_csv(C_file, sep='\t',header=None)
#l1_dist = compute_l1_distances(C)
#pred_labels = compute_k_means(l1_dist, n_clusters)
#save_clusters_file(pred_labels, out_file)

# getting true labels and clc ari
#true_labels = get_true_labels(Z_file)
#print(true_labels)
#kmeans_ari = adjusted_rand_score(true_labels, pred_labels)
#f = open(out_file+"/ARI_kmeans.txt",'w')
#f.write(str(kmeans_ari))
#f.close()


