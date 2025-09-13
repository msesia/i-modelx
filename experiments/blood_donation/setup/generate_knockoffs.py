import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import math
import csv
import itertools
import scipy.io
import pdb
from tqdm import tqdm
import pandas as pd
import sys
  
sys.path.append('metropolized_knockoffs')

seed = int(sys.argv[1])

N = 5

## Load the data
ifile = "../blood_donation_data/dt_donor_mar042018_simple.csv"
groups = np.loadtxt(ifile, delimiter="\t", skiprows=1, usecols=[0])

def group_to_treatments(groups):
    n = len(groups)
    X = -np.ones((n,N))
    
    idx_0 = np.where(groups==0)[0]
    X[idx_0,:] = [0,0,0,0,0]
    
    idx_1 = np.where(groups==1)[0]
    X[idx_1,:] = [1,0,0,0,0]
    
    idx_2 = np.where(groups==2)[0]
    X[idx_2,:] = [1,1,0,0,0]
    
    idx_3 = np.where(groups==3)[0]
    X[idx_3,:] = [1,0,1,0,0]
    
    idx_4 = np.where(groups==4)[0]
    X[idx_4,:] = [1,1,1,0,0]
    
    idx_5 = np.where(groups==5)[0]
    X[idx_5,:] = [1,0,1,1,0]

    idx_6 = np.where(groups==6)[0]
    X[idx_6,:] = [1,1,1,0,1]

    return X

X = group_to_treatments(groups)


# Log density
def lf(x):
    if (x[0]==0) and (x[1]==0) and (x[2]==0) and (x[3]==0) and (x[4]==0):      # T0 (no text)
        ans = 14/80
    elif (x[0]==1) and (x[1]==0) and (x[2]==0) and (x[3]==0) and (x[4]==0):    # T1 (text)
        ans = 11/80
    elif (x[0]==1) and (x[1]==1) and (x[2]==0) and (x[3]==0) and (x[4]==0):    # T2 (reward)
        ans = 11/80
    elif (x[0]==1) and (x[1]==0) and (x[2]==1) and (x[3]==0) and (x[4]==0):    # T3 (friend)
        ans = 11/80
    elif (x[0]==1) and (x[1]==1) and (x[2]==1) and (x[3]==0) and (x[4]==0):    # T4 (reward + friend)
        ans = 11/80
    elif (x[0]==1) and (x[1]==0) and (x[2]==1) and (x[3]==1) and (x[4]==0):    # T5 (friend + group)
        ans = 11/80
    elif (x[0]==1) and (x[1]==1) and (x[2]==1) and (x[3]==0) and (x[4]==1):    # T6 (reward + friend + small group)
        ans = 11/80
    else:
        ans = 1e-300
    return(ans)

def log_lf(x):
    return np.log(lf(x))

#sample from the given density function
def gibbs_sampler(lf, n_iter = 10000):
    d = N
    p_orig = 0
    # Make sure you start from a feasible configuration
    while p_orig < 1e-6:
        x_init = np.random.choice(2, size = d)
        p_orig = lf(x_init)
        
    x = x_init.copy()
    p_orig = lf(x)
    lks = [p_orig] #log probabilities, useful for diagnostics
    
    for _ in range(n_iter): 
        j = np.random.choice(d)
        orig = x[j]
        
        #metropolis update with uniformly random proposal
        prop = np.random.choice(2, size = 1)
        if prop == orig:
            continue      
        x[j] = prop
        p_prop = lf(x)    
        accept_prob = min(1, p_prop / np.maximum(1e-6,p_orig))
        u = np.random.uniform()
        if(u < accept_prob):
            p_orig = p_prop
            lks += [p_prop]
        else:
            x[j] = orig
            lks += [p_orig]
    
    return [x, x_init, lks]

#create adjacency matrix
adj_mat = np.ones([N,N])

#define the graph
G = nx.Graph()
G.add_nodes_from(range(N))
for i in range(N-1):
    for j in range(i+1,N):
        G.add_edge(i,j)


#tools to carry out the tree decomposition
import treewidth 
#tools to convert the tree decomposition to a format for metro
import graph_processing

def permute_order(order, frontier):
    new_order = np.random.permutation(order).tolist()
    N = len(order)
    new_active_frontier = []
    old_frontier = new_order
    for i in range(N):
        old_frontier = old_frontier[1:len(old_frontier)]
        new_active_frontier.append(old_frontier)
    return new_order, new_active_frontier

width, T = treewidth.treewidth_decomp(G) #get tree decomposition
#information for the metro sampler
order, active_frontier = graph_processing.get_ordering(T) 

import metro_generic as metro

#propose uniformly across states
def sym_proposal(j, xj):
    return np.random.choice(2)


## Generate the knockoffs
np.random.seed(seed)

n = X.shape[0]
nk = n
Xk = -np.ones((nk,N))
for i in tqdm(range(nk)):
    order, active_frontier = permute_order(order, active_frontier)
    Xk[i] = metro.single_metro(log_lf, X[i], order, active_frontier, sym_proposal)


def convert(x):
    try:
        return x.astype(int)
    except:
        return x

## Save the results
ofile = "../blood_donation_data/knockoffs/dt_donor_mar042018_simple_N"+str(N)+"_seed"+str(seed)+".csv"
df = pd.DataFrame({"X_1" : X[:,0], "X_2" : X[:,1], "X_3" : X[:,2], "X_4" : X[:,3], "X_5" : X[:,4],
                   "Xk_1" : Xk[:,0], "Xk_2" : Xk[:,1], "Xk_3" : Xk[:,2], "Xk_4" : Xk[:,3], "Xk_5" : Xk[:,4]})
df.apply(convert).to_csv(ofile, index=False, sep='\t', float_format=":d")

print("Done. Written knockoffs on {:s}".format(ofile))
