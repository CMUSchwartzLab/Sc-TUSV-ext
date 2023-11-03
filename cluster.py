import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib.pyplot as plt
import os
from scipy.spatial.distance import squareform
import sys


N_clusters = 2 * int(sys.argv[3]) -1  #  N_clusters: number of clones in the tree
filenames = [file.split('.')[0] for file in sorted(os.listdir(sys.argv[2]))]

arr_df = pd.read_csv(sys.argv[1]+'/input_pairwise_distances.tsv', sep='\t')

if len(filenames) >= 1 and 'diploid'<filenames[0]:
    arr = arr_df.values[1:,2:]#np.load(sys.argv[1]+'/input_pairwise_distances.tsv')[:-1,:-1]#('../medicc2_signatures_input/medicc2_output_DG1134_var_3/pairwise_distance.npy')[:-1,:-1]
else:
    arr = arr_df.values[:-1.1:-1]

dist_df = pd.DataFrame(arr,columns=filenames)
dist_df.set_index(filenames, drop=True, append=False, inplace=False, verify_integrity=False)
dist_df['sample_id']=filenames


dist_df.rename(columns = {'sample_id':'Cells'}, inplace = True)
cell_names = list(dist_df['Cells'])
distance_matrix = dist_df.iloc[:, :-1].values

# Normalize the distance matrix
normalized_distance_matrix = (distance_matrix - np.min(distance_matrix)) / (np.max(distance_matrix) - np.min(distance_matrix))
condensed_normalized_distance = normalized_distance_matrix[np.triu_indices(normalized_distance_matrix.shape[0], k=1)]

# Initial Hierarchical Clustering using "complete" linkage
Z_initial = linkage(condensed_normalized_distance, method='complete')
clusters_initial = fcluster(Z_initial, t=N_clusters, criterion='maxclust')



new_df = pd.DataFrame()
new_df['Cells'] = filenames
new_df['cluster'] = clusters_initial
new_df.to_csv(sys.argv[1]+'/clusters.tsv',sep='\t',index=False)
# Create a mapping from index to cluster label
index_to_cluster = {idx: cluster for idx, cluster in enumerate(clusters_initial)}

# Define custom leaf label function to show number of points in each leaf
def leaf_label_func(val):
    cluster = index_to_cluster[val]
    return '(%d points)' % (clusters_initial == cluster).sum()

# Plotting the dendrogram
fig, ax = plt.subplots()
#plt.figure(figsize=(20, 10))
d = dendrogram(Z_initial, 
           truncate_mode='lastp', 
           p=5, 
           show_leaf_counts=True, 
           show_contracted=True,
           color_threshold=5,
           leaf_font_size=14,
           orientation='right',
           ax = ax)

# Adjust the linewidths
for i, line in enumerate(ax.get_children()):
    if isinstance(line, plt.Line2D):
        line.set_linewidth(i/2 + 1)

plt.title("Hierarchical Clustering Dendrogram", fontsize=12)
plt.ylabel("Number of Points in Cluster", fontsize=12)
plt.xlabel("Distance", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.axhline(y=5, color='r', linestyle='--')  # Horizontal line indicating the cut-off for 5 clusters
plt.savefig(sys.argv[1]+'/cluster.png', dpi=300, bbox_inches='tight')
#plt.show()
