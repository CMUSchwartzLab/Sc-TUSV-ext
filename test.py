import copy
import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import random
import numpy as np
import multiprocessing as mp
import pickle
from datetime import datetime
from graphviz import Digraph
from ete2 import Tree          # for creating phylogenetic trees for .xml output
from Bio import Phylo          # for creating phylogenies to export as phylo .xml files
from cStringIO import StringIO # for converting string to file (for creating initial phylo .xml)
import pandas as pd
sys.path.insert(0, 'model/')
sys.path.insert(0, 'help/')
import solver as sv
import file_manager as fm      # sanitizes file and directory arguments
import generate_matrices as gm # gets F, Q, G, A, H from .vcf files
import printer as pt
import vcf_help as vh
import pickle
from snv_matching import snv_assign_sc_2
import cluster_l1 as cl



class ModifyTree:
    def __init__(self, E):
        self.cp_tree = {}
        self.tree = {}
        self.E = E
        self.N = len(E)
        for i in xrange(self.N - 1, -1, -1):
            for j in xrange(self.N - 1, -1, -1):
                if int(E[i, j]) == 1:
                    self.cp_tree[j] = i
                    if i not in self.tree.keys():
                        self.tree[i] = [j]
                    else:
                        self.tree[i].append(j)

    def delete_node(self, idx, Z):
        if self.is_root(idx):
            if self.num_children(idx) > 1:
                raise('Cannot delete root node with more than one child!')
            child = self.tree[idx][0]
            Z[:, child] = np.logical_or(Z[:, child], Z[:, idx]).astype(int) # NISHAT added: adding root's cells to its only child.
            Z[:, idx] = 0
            del self.cp_tree[child]
            del self.tree[idx]
        elif self.is_leaf(idx):
            parent = self.cp_tree[idx]
            del self.cp_tree[idx]
            if self.num_children(parent) == 1:
                del self.tree[parent]
                parents_parent = self.cp_tree[parent] # NISHAT added: if idx is leaf and idx is only child of its parent, map idx to its parent's parent.
                Z[:, parent] = np.logical_or(Z[:, parent], Z[:, idx]).astype(int) # NISHAT added: if idx is leaf and idx is only child of its parent, map idx to its parent's parent.
                Z[:, idx] = 0
            else:
                self.tree[parent].remove(idx)
                Z[:, parent] = np.logical_or(Z[:, parent], Z[:, idx]).astype(int) # NISHAT added: if idx is leaf and idx is not the only child of its parent, map idx to its parent.
                Z[:, idx] = 0
        else:
            parent = self.cp_tree[idx]
            children = self.tree[idx]
            del self.cp_tree[idx]
            for child in children:
                self.cp_tree[child] = parent
                self.tree[parent].append(child)
                self.E[parent, child] = 1
            del self.tree[idx]
            Z[:, parent] = np.logical_or(Z[:, parent], Z[:, idx]).astype(int) # NISHAT added: if idx is an internal node, map cells of idx to its parent.
            Z[:, idx] = 0
        return Z

    def is_leaf(self, idx):
        if idx not in self.tree.keys():
            return True
        else:
            return False

    def is_root(self, idx):
        if idx in self.tree.keys() and idx not in self.cp_tree.keys():
            return True
        else:
            return False

    def num_children(self, idx):
        if self.is_leaf(idx):
            return 0
        else:
            return len(self.tree[idx])
        
def collapse_nodes(C, E, A, R, Z, W, W_SV, W_SNV, threshold=0.0, only_leaf=False): # NISHAT: REMOVED U
    print("inside collapse")
    data_dict = {
    'C': pd.DataFrame(C),
    'E': pd.DataFrame(E),
    'A': pd.DataFrame(A),
    'R': pd.DataFrame(R),
    'Z': pd.DataFrame(Z),
    'W': pd.DataFrame(W),
    'W_SV': pd.DataFrame(W_SV),
    'W_SNV': pd.DataFrame(W_SNV),
    }

    # Example access: View matrix 'C'
    print("R", data_dict['R'])

    # Save the dictionary using pickle
    
    with open('matrices.pkl', 'wb') as f:
        pickle.dump(data_dict, f)
    # generate the tree
    tree = ModifyTree(E)
    if not only_leaf:
        # collapse the branches with 0 length
        branch_remove_idx = []
        for i in xrange(tree.N-1, -1, -1):
            for j in xrange(tree.N-1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and R[i,j] == 0:  # NISHAT: Added 'or sum(Z[:,j])==0' with the condition. 
                    branch_remove_idx.append(j)
        for j in xrange(tree.N-1, -1, -1):
            if sum(Z[:,j]) == 0:  # NISHAT: Added 'or sum(Z[:,j])==0' with the condition. 
                branch_remove_idx.append(j)
        branch_remove_idx = list(set(branch_remove_idx))
        print(tree.cp_tree)
        for node in branch_remove_idx:
            if node in tree.cp_tree.keys(): # if node to del is anythong except root
                print(node)
                target = tree.cp_tree[node] # get parent of node
                #U[:, target] += U[:, node]  
                if not tree.is_leaf(node): # if node not leaf
                    for child in tree.tree[node]: # for every child of node
                        print(child)
                        R[target, child] = R[node, child] # add node child's cost to parent's cost
                Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
                print("first",Z)

        # collapse the nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        '''  # NISHAT: REMOVED BLOCK
        for i in xrange(tree.N-1, -1, -1):
            if i in branch_remove_idx:
                continue
            if np.mean(U[:, i]) <= threshold:
                if tree.num_children(i) == 1:
                    freq_remove_idx.append(i)
                elif tree.is_leaf(i):
                    freq_leaf_remove_idx.append(i)
        '''
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
            print("2",Z)
            W[target, :] += W[node, :]
            W_SV[target,:] += W_SV[node,:]
            W_SNV[target, :] += W_SNV[node, :]
            R[parent, target] = R[parent, node] + R[node, target]
        for node in freq_leaf_remove_idx:
            Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
            print("3",Z)
    else:
        # collapse the branches with 0 length and the child of the branch doesn't belong to leaf nodes
        branch_remove_idx = []
        for i in xrange(tree.N - 1, -1, -1):
            for j in xrange(tree.N - 1, -1, -1):
                if int(E[i, j]) == 1 and sum(W[j, :]) == 0 and R[i, j] == 0 and not tree.is_leaf(j):
                    branch_remove_idx.append(j)
        for node in branch_remove_idx:
            target = tree.cp_tree[node]
            #U[:, target] += U[:, node] # NISHAT: REMOVED
            if not tree.is_leaf(node):
                for child in tree.tree[node]:
                    R[target, child] = R[node, child]
            Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
            print("4",Z)

        # collapse the leaf nodes with 0 frequency
        freq_remove_idx = []
        freq_leaf_remove_idx = []
        for i in xrange(tree.N - 1, -1, -1):
            if i in branch_remove_idx:
                continue
            #if np.mean(U[:, i]) <= threshold and tree.is_leaf(i): # NISHAT: REMOVED U's condition
            if tree.is_leaf(i):
                freq_leaf_remove_idx.append(i)
        for node in freq_leaf_remove_idx:
            Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
            print(Z)
        for i in xrange(tree.N - 1, -1, -1):
            if tree.num_children(i) == 1:
                freq_remove_idx.append(i)
        for node in freq_remove_idx:
            target = tree.tree[node][0]
            parent = tree.cp_tree[node]
            Z = tree.delete_node(node, Z) # NISHAT: Added Z to this line. Previous: tree.delete_node(node)
            print("5",Z)
            W[target, :] += W[node, :]
            W_SV[target, :] += W_SV[node, :]
            W_SNV[target, :] += W_SNV[node, :]
            R[parent, target] = R[parent, node] + R[node, target]


    # delete those nodes
    remove_idx = branch_remove_idx + freq_remove_idx + freq_leaf_remove_idx
    print('Nodes ', remove_idx, 'will be collapsed.')
    #U_new = np.delete(U, remove_idx, axis=1) # NISHAT
    C_new = np.delete(C, remove_idx, axis=0)
    A_new = np.delete(A, remove_idx, axis=0)
    A_new = np.delete(A_new, remove_idx, axis=1)
    E_new = np.delete(tree.E, remove_idx, axis=0)
    E_new = np.delete(E_new, remove_idx, axis=1)
    R_new = np.delete(R, remove_idx, axis=0)
    R_new = np.delete(R_new, remove_idx, axis=1)
    Z_new = np.delete(Z, remove_idx, axis=1) # NISHAT: added for Z.
    print(Z_new)
    W_new = np.delete(W, remove_idx, axis=0)
    W_SV_new = np.delete(W_SV, remove_idx, axis=0)
    W_SNV_new = np.delete(W_SNV, remove_idx, axis=0)
    #print("collapse", U_new.shape, C_new.shape) # NISHAT: REMOVED U
    # return U_new, C_new, E_new, A_new, R_new, W_new, W_SV_new, W_SNV_new # NISHAT: REMOVED U
    return C_new, E_new, A_new, R_new, Z_new, W_new, W_SV_new, W_SNV_new, remove_idx


with open('matrices.pkl', 'rb') as f:
    loaded_data = pickle.load(f)

# Access a matrix
C = loaded_data['C'].values
E = loaded_data['E'].values
A = loaded_data['A'].values
R = loaded_data['R'].values
Z = loaded_data['Z'].values
W = loaded_data['W'].values
W_SNV = loaded_data['W_SNV'].values
W_SV = loaded_data['W_SV'].values
print(Z)
C_new, E_new, A_new, R_new, Z_new, W_new, W_SV_new, W_SNV_new, remove_idx = collapse_nodes(C, E, A, R, Z, W, W_SV, W_SNV, threshold=0.0, only_leaf=False)
print(Z_new)