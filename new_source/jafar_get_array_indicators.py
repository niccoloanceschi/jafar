import numpy as np

def make_tensor_py(M, K):
    M = int(M)
    K = int(K)
    shape = [K] * (M + 1)
    grids = np.meshgrid(*[np.arange(1, K + 1)] * (M + 1), indexing='ij')

    i1 = grids[0]
    h = grids[-1]
    rest = grids[1:M]

    max_rest = np.maximum.reduce(rest)
    mask = (i1 > h) & (max_rest > h)
    return mask.astype(np.int8)




