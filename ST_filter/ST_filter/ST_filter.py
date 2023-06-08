"""
Identifying significant ties by fitness model (Kobayashi et.al, 2019).

Reference: 
T.Kobayashi, T.Takaguchi, A.Barrat (2019) "The structured backbone of temporal 
    social ties", Nature Communications 10.
----------------------------------------------------------------------------
This program is free software following MIT License.

Copyright (c) 2023 Yoshitaka Ogisu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights 
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.

yoshitaka.ogisu@gmail.com
https://sites.google.com/view/yoshitakaogisu/
----------------------------------------------------------------------------
"""


import numpy as np
import scipy.stats as stat
import scipy.optimize as opt
from scipy.sparse import csr_matrix
from itertools import combinations
import time

def ST_filter(A_snap, alpha, 
              judge='p_val', opt_method='krylov'):
    """
    Identifying significant ties from list of adjacency matrix.
    
    Input 
    -----
    A_snap : 3D-array_like (t, i, j)
        Snapshots of undirected matrix at time t.
        Each matrix contains only 1 or 0 for every element.
    alpha : float
        Significant level.
    judge : str (default is 'p_val')
        What significant ties are judged based on.
        'p_val' : p-values calculated by activity parameters
        'inv_binom' : inverse function of binomial function
    opt_method : str (default is 'krylov')
        Optimization method for nonlinear root finding probrem.
        See method options in 'scipy.optimize.root'.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
    
    Output
    ------
    Adj_sig : 2D-array_like
        Matrix in which significant and insignificant ties are 1 and 0, respectively.
    Adj_all : 2D-array_like
        Adjacency matrix of the network summing up overall periods.
    p_mat : 2D-array_like
        Matrix showing the p-value for each edge in statistical test.
    activ_params : array_like
        Estimated activity parameter of each node.    
    """
    
    tau = len(A_snap) # number of snapshots
    N = len(A_snap[0]) # number of nodes for each matrix
    
    Adj_all = np.sum(A_snap,axis=0) # Adjacency matrix summing up for all time
    #Adj_all = csr_matrix(np.sum(A_snap,axis=0)) # sparse ver. Adjacency matrix summing up for all time
    
    def obj_func(x):
        """
        Function for optimization
        """
        x = x.reshape(-1,1)
        ## probability function
        u = x@x.T
        obj_mat = (Adj_all - tau*u)/(1-u)
        H = np.sum(obj_mat, axis=1) - np.diag(obj_mat) # sum up i neq j
        
        return H
    
    ## configuration model for first guess
    numer = np.sum(Adj_all, axis=1)/tau
    denom = np.sqrt(np.sum(Adj_all)/tau)
    a0 = numer/denom
    
    ## estimate activity parameters
    sol = opt.root(obj_func, x0=a0, method=opt_method)
    if sol.success == True:
        print("Root finding method completed successfully.")
    else:
        print("Caution! Root finding method failed to find an appropriate solution.")
        print("The solution could be wrong.")
    activ_params = sol.x
    
    ## construct probability matrix
    u = activ_params.reshape(-1,1)@activ_params.reshape(1,-1)
    
    #Adj_all = Adj_all.toarray()
    ## calculate p-values
    p_mat = stat.binom.sf(k = Adj_all, n = tau, p = u)+np.diag(np.nan+np.ones(len(Adj_all))) # diagonal is nan
    
    ## compute significant edges
    if judge == 'p_val':
        sigs = np.array(p_mat <= alpha, int)
        Adj_sig = np.array(Adj_all>0,int)*sigs # set 0 if each node pair has no edge in Adj_all
    elif judge == 'inv_binom':
        Adj_sig = np.array(Adj_all > stat.binom.isf(q = alpha, n=tau, p=u), int)
    else:
        Adj_sig = np.array(p_mat <= alpha, int)
        print("Caution! 'judge' can take only 'p_val' or 'inv_binom'.")
        print("Now we calculate significant ties based on p-values.")
        
        
    
    return {"Adj_sig": Adj_sig,
            "Adj_all": Adj_all, 
            "p_mat": p_mat, 
            "activ_params": activ_params}


def ST_filter_aat(AAT, t, alpha, 
              judge='p_val', opt_method='krylov'):
    """
    Identifying significant ties from aggregate adjacency matrix.
    
    Input 
    -----
    AAT : 2D-array_like (i, j)
        Aggregate adjacency matrix (undirected).
        Each matrix contains only 0 or nonzero integer for every element.
    t: number of snapshots.
    alpha : float
        Significant level.
    judge : str (default is 'p_val')
        What significant ties are judged based on.
        'p_val' : p-values calculated by activity parameters
        'inv_binom' : inverse function of binomial function
    opt_method : str (default is 'krylov')
        Optimization method for nonlinear root finding probrem.
        See method options in 'scipy.optimize.root'.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
    
    Output
    ------
    Adj_sig : 2D-array_like
        Matrix in which significant and insignificant ties are 1 and 0, respectively.
    Adj_all : 2D-array_like
        Adjacency matrix of the network summing up overall periods.
    p_mat : 2D-array_like
        Matrix showing the p-value for each edge in statistical test.
    activ_params : array_like
        Estimated activity parameter of each node.    
    """
    
    tau = t # number of snapshots
    N = len(AAT) # number of nodes for each matrix
    
    Adj_all = AAT # Adjacency matrix summing up for all time
    #Adj_all = csr_matrix(np.sum(A_snap,axis=0)) # sparse ver. Adjacency matrix summing up for all time
    
    def obj_func(x):
        """
        Function for optimization
        """
        x = x.reshape(-1,1)
        ## probability function
        u = x@x.T
        obj_mat = (Adj_all - tau*u)/(1-u)
        H = np.sum(obj_mat, axis=1) - np.diag(obj_mat) # sum up i neq j
        
        return H
    
    ## configuration model for first guess
    numer = np.sum(Adj_all, axis=1)/tau
    denom = np.sqrt(np.sum(Adj_all)/tau)
    a0 = numer/denom
    
    ## estimate activity parameters
    sol = opt.root(obj_func, x0=a0, method=opt_method)
    if sol.success == True:
        print("Root finding method completed successfully.")
    else:
        print("Caution! Root finding method failed to find an appropriate solution.")
        print("The solution could be wrong.")
    activ_params = sol.x
    
    ## construct probability matrix
    u = activ_params.reshape(-1,1)@activ_params.reshape(1,-1)
    
    #Adj_all = Adj_all.toarray()
    ## calculate p-values
    p_mat = stat.binom.sf(k = Adj_all, n = tau, p = u)+np.diag(np.nan+np.ones(len(Adj_all))) # diagonal is nan
    
    ## compute significant edges
    if judge == 'p_val':
        sigs = np.array(p_mat < alpha, int)
        Adj_sig = np.array(Adj_all>0,int)*sigs # set 0 if each node pair has no edge in Adj_all
    elif judge == 'inv_binom':
        Adj_sig = np.array(Adj_all > stat.binom.isf(q = alpha, n=tau, p=u), int)
    else:
        Adj_sig = np.array(p_mat < alpha, int)
        print("Caution! 'judge' can take only 'p_val' or 'inv_binom'.")
        print("Now we calculate significant ties based on p-values.")
        
        
    
    return {"Adj_sig": Adj_sig,
            "Adj_all": Adj_all, 
            "p_mat": p_mat, 
            "activ_params": activ_params}

def ST_filter_list(edge_list, alpha, 
                   judge='p_val', opt_method='krylov',
                   memorysave=False, paral=True):
    """
    Identifying significant ties from edge list.
    
    Input 
    -----
    edge_list : 2D-array_like
        Edge list of temporal network.
        Each row of the list needs to organize such that ["Snapshot ID", "node1", "node2"].
    alpha : float
        Significant level.
    judge : str (default is 'p_val')
        What significant ties are judged based on.
        'p_val' : p-values calculated by activity parameters
        'inv_binom' : inverse function of binomial function
    opt_method : str (default is 'krylov')
        Optimization method for nonlinear root finding probrem.
        See method options in 'scipy.optimize.root'.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html#scipy.optimize.root
    memorysave : bool (default is 'False')
        If memory error arised, you can choose memory save mode.
        Note1 : if 'memorysave == True', the algorithm is much slower.
        Note2 : When 'memorysave == True', the return value is 'activ_params' only.
    paral : bool (default is 'True')
        When using memory save mode, parallel computing by threading process is used.
    
    
    Output
    ------
    Adj_sig : 2D-array_like
        Matrix in which significant and insignificant ties are 1 and 0, respectively.
    Adj_all : 2D-array_like
        Adjacency matrix of the network summing up overall periods.
    p_mat : 2D-array_like
        Matrix showing the p-value for each edge in statistical test.
    activ_params : array_like
        Estimated activity parameter of each node.
    
    Output (Additional)
    -------------------
    nodes : array_like (only when edge_list includes str objects)
        Nodes list of the network
        
    """
    edge_list = np.array(edge_list)
    
    tau = len(set(edge_list.T[0])) # number of snapshots
    nodes = list(set(list(set(edge_list.T[1]))+list(set(edge_list.T[2]))))
    N = len(nodes) # number of nodes for each matrix
    
    # allocate node id if including str
    try:
        edge_list = np.array(edge_list, int)
        in_str = False
    except ValueError: # if edge_list includes strings, replace them to numbers
        for i in range(N):
            edge_list = np.where(edge_list==nodes[i], i, edge_list)
        edge_list = np.array(edge_list,int) # transrate dtype to integer
        in_str = True
        
        
    if memorysave==False:
        try:        
            ## try algorithm being fast but huge memory needed
            e_list = edge_list[:,[1,2]]
            data= np.ones(len(e_list))
            u_agg_mat = csr_matrix((data, (e_list.T)),shape=(N,N))
            Adj_all = (u_agg_mat+u_agg_mat.T).toarray()

            def obj_func(x):
                """
                Function for optimization
                """
                x = x.reshape(-1,1)
                ## probability function
                u = x@x.T
                obj_mat = (Adj_all - tau*u)/(1-u)
                H = np.sum(obj_mat, axis=1) - np.diag(obj_mat) # sum up i neq j

                return H

            ## configuration model for first guess
            numer = np.sum(Adj_all, axis=1)/tau
            denom = np.sqrt(np.sum(Adj_all)/tau)
            a0 = numer/denom

            ## estimate activity parameters
            sol = opt.root(obj_func, x0=a0, method=opt_method)
            if sol.success == True:
                print("Root finding method completed successfully.")
            else:
                print("Caution! Root finding method failed to find an appropriate solution.")
                print("The solution could be wrong.")
            activ_params = sol.x
            ## construct probability matrix
            u = activ_params.reshape(-1,1)@activ_params.reshape(1,-1)

            ## calculate p-values
            p_mat = stat.binom.sf(k = Adj_all, n = tau, p = u)+np.diag(np.nan+np.ones(len(Adj_all))) # diagonal is nan

            ## compute significant edges
            if judge == 'p_val':
                sigs = np.array(p_mat <= alpha, int)
                Adj_sig = np.array(Adj_all>0,int)*sigs # set 0 if each node pair has no edge in Adj_all
            elif judge == 'inv_binom':
                Adj_sig = np.array(Adj_all > stat.binom.isf(q = alpha, n=tau, p=u), int)
            else:
                Adj_sig = np.array(p_mat < alpha, int)
                print("Caution! 'judge' can take only 'p_val' or 'inv_binom'.")
                print("Now we calculate significant ties based on p-values.")
            
            if in_str == False:
                resultset = {"Adj_sig": Adj_sig,
                             "Adj_all": Adj_all, 
                             "p_mat": p_mat, 
                             "activ_params": activ_params}
            else:
                resultset = {"Adj_sig": Adj_sig,
                             "Adj_all": Adj_all, 
                             "p_mat": p_mat, 
                             "activ_params": activ_params,
                             "nodes": nodes}
    
        except MemoryError:
            print('"MemoryError" arised.')
            print('Try memory saving algorithm "memorysave=True" if you wish.')
            # end the algorithm
            return
        
    else: # memory save mode
        print('Caution! Now using memory saving mode.')
        print('The run time is much slower than the default algorithm.')
        ## if MemoryError appears, use memory saving algolithm 
        search_list = np.sort(edge_list[:,[1,2]],axis=1)

        if paral == True:
            from concurrent.futures import ThreadPoolExecutor
            def obj_func(x):
                """
                Function for optimization
                """ 
                H = np.zeros(N)
                def sol(i):
                    ## probability function
                    u = x[i]*x
                    edges = search_list[np.where(search_list==i)[0]]
                    connected_nodes = np.array(list(set(edges.reshape(-1))))
                    mij = np.histogram(edges[np.where(edges!=i)],
                                 bins=N,
                                 range=(0,N) )[0]
                    res = (mij-tau*u)/(1-u)
                    res[i] = 0
                    return res.sum()
                with ThreadPoolExecutor() as executor:
                    result = executor.map(sol,range(N))
                    H = np.array(list(result))

                return H
            ## configuration model for first guess
            numer = np.array([np.sum(edge_list[:,[1,2]]==j) for j in range(len(nodes))])/tau
            denom = np.sqrt(len(edge_list)*2/tau)
            a0 = numer/denom


            ## estimate activity parameters
            sol = opt.root(obj_func, x0=a0, method=opt_method)
            if sol.success == True:
                print("Root finding method completed successfully.")
            else:
                print("Caution! Root finding method failed to find an appropriate solution.")
                print("The solution could be wrong.")
                
        else:
            def obj_func(x):
                """
                Function for optimization
                """
                H = np.zeros(N)
                for i in range(N):
                    ## probability function
                    u = x[i]*x
                    edges = search_list[np.where(search_list==i)[0]]
                    connected_nodes = np.array(list(set(edges.reshape(-1))))
                    mij = np.histogram(edges[np.where(edges!=i)],
                                 bins=N,
                                 range=(0,N) )[0]
                    res = (mij-tau*u)/(1-u)
                    res[i] = 0
                    H[i] = res.sum()

                return H    
            ## configuration model for first guess
            numer = np.array([np.sum(edge_list[:,[1,2]]==j) for j in range(len(nodes))])/tau
            denom = np.sqrt(len(edge_list)*2/tau)
            a0 = numer/denom


            ## estimate activity parameters
            sol = opt.root(obj_func, x0=a0, method=opt_method)
            if sol.success == True:
                print("Root finding method completed successfully.")
            else:
                print("Caution! Root finding method failed to find an appropriate solution.")
                print("The solution could be wrong.")
        
        activ_params = sol.x
        if in_str==False:
            resultset = {"activ_params": activ_params}
        else:
            resultset = {"activ_params": activ_params,
                         "nodes": nodes}
        
    return resultset

