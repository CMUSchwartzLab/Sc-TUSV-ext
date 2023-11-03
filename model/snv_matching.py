#  author: Xuecong Fu
#  Match unsampled SNVs (or even breakpoints if necessary) to inferred phylogeny given the results from TUSV-ext

import numpy as np
import re


def dot2pctable(dotfile):
    parent_child_table = []
    file = open(dotfile, encoding='utf8')
    line = file.readline().strip()
    if len(line) > 1:
        line = file.readline().strip()
        match = re.search(r'\d+\s->\s\d+.*', line)
        if match != None:
            pc_pair = re.match(r'(\d+)\s->\s(\d+).*(\d)\.',line).groups()
            parent_child_table.append(np.array(list(pc_pair)).astype(int))
    return np.array(parent_child_table)

def one_hot_encoding(index_list, category_list):
    encoding = np.zeros((len(category_list), len(index_list)))
    encoding[:, index_list] = 1
    return encoding

def get_diff(E,C_SNV_avg):
    
    C_SNV_avg = C_SNV_avg.reshape((E.shape[0],1))
    dist = np.abs(C_SNV_avg*E - C_SNV_avg.transpose()*E)
    
    i,j = np.unravel_index(np.argmax(dist), dist.shape)
    return j, np.max(dist)


def snv_assign_sc(C_SNV_unsampled, Z, E):
    #print('mapping unsampled snvs')
    N, n = Z.shape
    if C_SNV_unsampled != None:
        l_g_un = C_SNV_unsampled.shape[1]
    else:
        l_g_un = 0
    
    min_dist = np.full((l_g_un), 10000)
    min_node = np.full((l_g_un), -1)
    
    W_snv = np.zeros((n,l_g_un))
    
    indices = [np.where(Z[:, i] == 1)[0] for i in range(n)]
    C_SNV_avg = np.array([np.sum(C_SNV_unsampled[indices[i]], axis=0)/len(indices[i]) for i in range(n)])
    
    C_SNV_avg[np.isnan(C_SNV_avg)] = 0
    #print(C_SNV_avg)
    
    for i in range(l_g_un): 
        min_node[i], min_dist[i] = get_diff(E,C_SNV_avg[:,i])
        W_snv[min_node[i], i] = 1
    return min_node, min_dist, W_snv

def snv_assign_sc_2(C_SNV_unsampled, Z, E, A):
    #print('Mapping unsampled snvs and SVs')
    """
    the function for assigning unsampled SNVs to the trees, with minimum
    distance criteria to identify the possible branch of a SNV and SV given
    n - number of clones
    N - number of cells
    l_g_un - number of unsampled variants
    :param C_SNV_unsampled: N*l_g_un unsampled SNV/SV for each cell (from input data)
    :param C_SNV_unsampled_clone: n*l_g_un unsampled SNV for each clone/node (calculated by majority voting)
    :param C_expected: n*l_g_un, c_snv_expected_i unsampled SNV for each clone if SNV/SV acquaired in node i.
    :param A: n*n, a_ij = 1 if i is the ancestor of j, diagonal is 0, which means i is not the ancestor of i
    :return: node and edge assignments of unsampled SNVs and SVs.
    """
    N, n = Z.shape
    if not np.isnan(C_SNV_unsampled).any():
        l_g_un = C_SNV_unsampled.shape[1]
    else:
        l_g_un = 0
    
    # Convert copy numbers of SNVs and SVs to 1s or 0s. 1s if copy number > 0.
    C_SNV_unsampled = (C_SNV_unsampled>0)*1
    
    # Assign clone's/node's SNV or SV according to majority voting of the cells in the clone/node.
    C_SNV_unsampled_clone = np.array([ # NISHAT COMMENTED oct: We do not need to shuffle C's rows according to Z as we are mapping to the tree. C_est is already in the order of the tree. 
        np.apply_along_axis(lambda x: np.argmax(np.bincount(x)) if len(x) > 0 else np.nan, axis=0, 
                            arr=C_SNV_unsampled[np.argwhere(Z[:, i] == 1)])[0] for i in range(n)
                            ])  
    #C_SNV_unsampled_clone = C_SNV_unsampled
    
    C_expected = A + np.eye(n)
    
    

    min_node = np.full((l_g_un), -1)
    dists = np.zeros((n,l_g_un))
    
    dists = [np.nansum(np.abs(C_expected[i,:].reshape(n,1) - C_SNV_unsampled_clone),axis=0) for i in range(n)]
    
    # Compute the minimum value for each snv/sv pos
    min_dist = np.amin(dists, axis=0)

    # Compute the index of the minimum value for each snv/sv pos
    min_node = np.argmin(dists, axis=0)
    
    
    W_snv = np.zeros((n,l_g_un))
    for i in range(l_g_un): 
        W_snv[min_node[i], i] = 1
    return min_node, min_dist, W_snv

def snv_assign(C_CNV, Q, A, E, U, F, G):
#def snv_assign(C_CNV, Q, A, E, F, G): # Nishat removed U
    """
    the function for assigning unsampled SNVs to the trees, using brutal force with minimum
    distance criteria to identify the possible branch and allele of a SNV given
    n - number of clones
    m - number of samples
    l - number of SVs
    l_un - number of unsampled SVs
    g - number of sampled SNVs
    g_un - number of unsampled SNVs
    r - number of CNVs
    :param C_CNV: n*2r allelic specific CNV
    :param Q: (l_un + g_un) * r mapping matrix which maps the unsampled SNVs to CNV segments, q_ij=1 if ith SNV maps to jth CNV
    :param A: n*n, a_ij = 1 if i is the ancestor of j, diagonal is 0, which means i is not the ancestor of i
    :param U: m*n frequency matrix
    :param F: m*g_un frequency matrix
    :param G: l_un*l_un unsampled breakpoints pairing matrix
    :return:
    """
    n, r = C_CNV.shape
    l_g_un = Q.shape[0]
    l_un = G.shape[0]
    r = int(r/2)
    clone_idx_range = range(0, n-1) 
     # exclude the root node
    C_hat_1 = np.dot(C_CNV[:, :r], np.transpose(Q)) # n*l_g_un, the copy number of CNV at SNV position
    C_hat_2 = np.dot(C_CNV[:, r:], np.transpose(Q))  # n*l_g_un, the copy number of CNV at SNV position
    C_hat_1_parent = np.dot(E.T, C_hat_1)
    C_hat_2_parent = np.dot(E.T, C_hat_2)
    min_dist = np.full((l_g_un), 10000)
    min_node = np.full((l_g_un), -1)
    dist = np.full((l_g_un), np.inf)
    for b in clone_idx_range:
        ### normal copy number=1
        C_SNV_clone_1 = C_hat_1[b, :] # l_g_un
        C_SNV_clone_2 = C_hat_2[b, :]
        C_SNV_clone_parent_1 = C_hat_1_parent[b, :]
        C_SNV_clone_parent_2 = C_hat_2_parent[b, :]

        valid_snv_idx = np.array(list(set(np.append(np.where(C_SNV_clone_1 == 1)[0],np.where(C_SNV_clone_1 - C_SNV_clone_parent_1 > 1)[0]))))
        F_est = U[:,b][:,np.newaxis] * C_SNV_clone_1[valid_snv_idx] + np.dot(U, A[b, :][:,np.newaxis]* C_hat_1[:,valid_snv_idx])
        dist[valid_snv_idx] = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)
        #dist_stack = np.column_stack((min_dist[valid_snv_idx], dist))
        # argmin = np.argmin(dist_stack, axis=-1)
        # if (argmin == 1).any():
        #     #print(min_node[valid_snv_idx[argmin == 1]])
        #     min_node[valid_snv_idx[argmin == 1]] = b
        #     #print(min_node[valid_snv_idx[argmin == 1]])
        #     min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
        #     #print(min_node, min_dist)

        ### copy number > 1
        valid_snv_idx2 = np.where(C_SNV_clone_1 > 1)[0]
        F_est = U[:, b][:, np.newaxis] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_1[:, valid_snv_idx2] / C_SNV_clone_1[
                                                    valid_snv_idx2])
        dist[valid_snv_idx2] = np.sum(np.abs(F_est - F[:, valid_snv_idx2]), axis=0)

        dist[: l_un] += np.dot(dist[:l_un], G) #add the other corresponding breakpoint distance to original breakpoint to ensure paired breakpoints are at the same node
        dist[: l_un] /= 2
        dist_stack = np.column_stack((min_dist, dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            min_node[argmin == 1] = b
            min_dist = np.min(dist_stack, axis=-1)

        valid_snv_idx = np.array(list(set(np.append(np.where(C_SNV_clone_2 == 1)[0],np.where(C_SNV_clone_2 - C_SNV_clone_parent_2 > 1)[0]))))
        F_est = U[:, b][:,np.newaxis] * C_SNV_clone_2[valid_snv_idx] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_2[:, valid_snv_idx])
        dist[valid_snv_idx] = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)

        valid_snv_idx2 = np.where(C_SNV_clone_2 > 1)[0]
        F_est = U[:, b][:,np.newaxis] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_2[:, valid_snv_idx2] / C_SNV_clone_2[valid_snv_idx2])
        dist[valid_snv_idx2] = np.sum(np.abs(F_est - F[:, valid_snv_idx2]),axis=0)

        dist[: l_un] += np.dot(dist[:l_un], G) #add the other corresponding breakpoint distance to original breakpoint to ensure paired breakpoints are at the same node
        dist[: l_un] /= 2

        dist_stack = np.column_stack((min_dist, dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            min_node[valid_snv_idx[argmin == 1]] = b
            min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
    W_snv = np.zeros((n, len(min_node)))
    for i in range(len(min_node)):
        W_snv[min_node[i], i] = 1
    return min_node, min_dist, W_snv



if __name__ == '__main__':
    ### test case
    n=3
    m=1
    g_un=5
    k=4
    #np.random.seed(0)
    U = np.array([[0.1, 0.5, 0.3, 0.1],[0.2, 0.2, 0.6, 0]])
    C_CNV = np.array([[1,2,1,1,1,1,4,1],
                      [2,2,1,3,1,1,1,1],
                      [1,2,1,1,1,1,1,1],
                      [1,1,1,1,1,1,1,1],])
    A = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [1, 1, 0, 0],
                  [1, 1, 1, 0]])
    Q = np.eye(4)
    G = np.eye(2)
    C_SNV = np.array([[1, 1, 4, 0],
                      [2, 1, 1, 1],
                      [1, 1, 1, 0],
                      [0, 0, 0, 0]])
    E = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [1, 1, 0, 0],
                  [0, 0, 1, 0]])
    F_true = np.dot(U, C_SNV)
    F_noise = F_true + np.random.normal(scale=0.2,size=np.shape(F_true) )
    min_node, min_dist, W_snv = snv_assign(C_CNV, Q, A, E, U, F_noise, G)
    
