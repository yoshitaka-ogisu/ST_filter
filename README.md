# ST filter for python

## 1. Abstract
Find *significant ties* on unweited and undirected temporal networks.

## 2. Environment (test environment)
* Windows 11, Mac OS 12
* Python 3.7, 3.8
*
## 3. Dependency
* Numpy
* Scipy

## 4. Usage
Install by local pip or put **ST_filter.py** on working directory, and
```python
import ST_filter
```

## 5. Functions
```python
def ST_filter(A_snap, alpha, judge='p_val', opt_method='krylov')
```

### inputs

`A_snap` A list of snapshots as adjacency matrices. Each sanpshot is $n \times n$.

`alpha` Significant level (e.g. 0.1, 0.05, 0.01).

`judge` (default is `p_val`) What is the basis for detecting significant ties?
* `p_val` p-values < alpha
* `inv_binom` number of edges based on inverse binomial function < observed number of edges

`opt_method` (default is `krylov`) Algorighm to solve nonlinear system of equations. See method in [scipy.optimize.root](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root.html).

### outputs

`Adj_sig` $n \times n$ matrix in which an element is 1 if it is a significant tie, otherwise 0.

`Adj_all` $n \times n$ matrix in which each element is the observed number of edges during the sample periods.

`p_mat` p-values for observed number of edges.

`activ_params` avtivity parameters for each nodes.


 
