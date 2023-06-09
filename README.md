# ST filter for python

## 1. Abstract
Find *significant ties* on unweited and undirected temporal networks.

## 2. Environment (test environment)
* Windows 11, Mac OS 12
* Python 3.7, 3.8

## 3. Dependency
* Numpy
* Scipy

## 4. Usage
Install by local pip or put **ST_filter.py** on working directory, and
```python
import ST_filter
```
See **7. Installation**

## 5. Function

```python
def ST_filter(A_snap, alpha, judge='p_val', opt_method='krylov')
```

**inputs**
>`A_snap` List of snapshots as adjacency matrices. Each sanpshot is $n \times n$.
>
>`alpha` Significant level (e.g. 0.1, 0.05, 0.01).
>
>`judge` (default is `p_val`) What is the basis for detecting significant ties?
>* `p_val` p-values < alpha
>* `inv_binom` number of edges based on inverse binomial function < observed number of edges
>
>`opt_method` (default is `krylov`) Algorighm to solve nonlinear system of equations. See method in [scipy.optimize.root](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html).

**outputs**
>`Adj_sig` $n \times n$ matrix in which an element is 1 if it is a significant tie, otherwise 0.
>
>`Adj_all` $n \times n$ matrix in which each element is the observed number of edges during the sample periods.
>
>`p_mat` p-values for observed number of edges.
>
>`activ_params` avtivity parameters for each nodes.

---

```python
def ST_filter_list(edge_list, alpha, judge='p_val', opt_method='krylov')
```
**inputs**
>`edge_list` Edge list of a temporal network. Each row of the list consists `['Snapshot_ID', 'Node_ID#1', 'Node_ID#2']`.
> `Node_ID` can be strings. 
> `[x, 'Node_ID#1', 'Node_ID#2']` and `[x, 'Node_ID#2', 'Node_ID#1']` is treated as the same one.
>
>`alpha` Significant level (e.g. 0.1, 0.05, 0.01).
>
>`judge` (default is `p_val`) What is the basis for detecting significant ties?
>* `p_val` p-values < alpha
>* `inv_binom` number of edges based on inverse binomial function < observed number of edges
>
>`opt_method` (default is `krylov`) Algorighm to solve nonlinear system of equations. See method in [scipy.optimize.root](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html).
 
**outputs**
>`Adj_sig` $n \times n$ matrix in which an element is 1 if it is a significant tie, otherwise 0.
>
>`Adj_all` $n \times n$ matrix in which each element is the observed number of edges during the sample periods.
>
>`p_mat` p-values for observed number of edges.
>
>`activ_params` avtivity parameters for each nodes.
>
>`nodes` (rertun if `Node_ID` is given strings) Name of nodes which is listed in order of correspondence to activ_params.

---

```python
def ST_filter_aat(AAT, t, alpha, judge='p_val', opt_method='krylov')
```
**inputs**
>`AAT` Aggregate adjacency matrix (undirected) which is sum of snapshots of adjacency matrices.
>
>`t` Number of snapshots.
>
>`alpha` Significant level (e.g. 0.1, 0.05, 0.01).
>
>`judge` (default is `p_val`) What is the basis for detecting significant ties?
>* `p_val` p-values < alpha
>* `inv_binom` number of edges based on inverse binomial function < observed number of edges
>
>`opt_method` (default is `krylov`) Algorighm to solve nonlinear system of equations. See method in [scipy.optimize.root](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html).
 
**outputs**
>`Adj_sig` $n \times n$ matrix in which an element is 1 if it is a significant tie, otherwise 0.
>
>`Adj_all` $n \times n$ matrix in which each element is the observed number of edges during the sample periods.
>
>`p_mat` p-values for observed number of edges.
>
>`activ_params` avtivity parameters for each nodes.

## 6.Example
*Example 1*
```python
import ST_filter as stf
# Make test matrix A having 10 snapshots and 20 nodes.
A = [np.random.choice(a = [0,1], size = 400, p= [0.9,0.1]).reshape((20,20)) for i in range(10)]
for t in range(len(A)):
    for i in range(len(A[t])):
        for j  in range(i,len(A[t])):
            if i == j:
                A[t][i][j] = 0
            else:
                A[t][j][i] = A[t][i][j]
                
result = stf.ST_filter(A, 0.1) # use ST_filter
print(result)
```

*Example 2*
```python
from ST_filter import *

# Make test list E0.
E0 = np.array([[0, 2, 3],
                [0, 5, 10],
                [1, 1, 2],
                [1, 3, 4],
                [2, 2, 10]]) # [SnapshotID, node1, node2] for each row

result = ST_filter_list(E0, 0.1) # use ST_filter_list
print(result)
```

## 7. Installation
