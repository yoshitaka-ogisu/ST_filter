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
### input
A_snap
A list of snapshots as adjacency matrices. Each sanpshot is $n \times n$.
 
