# ===============================================================================================================
# SOURCE: https://github.com/pbansal5/DeepMVI
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/abs/2103.01600
# ===============================================================================================================

import numpy as np
import copy
from contextlib import contextmanager

@contextmanager
def null_context():
    yield

def is_blackout(matrix):
    arr = (np.sum(np.isnan(matrix).astype(int), axis=1) == matrix.shape[1])
    return arr.astype(int).sum() > 0




def make_validation (matrix,num_missing=20):
    np.random.seed(0)
    nan_mask = np.isnan(matrix)

    padded_mat = np.concatenate([np.zeros((1,nan_mask.shape[1])), nan_mask,np.zeros((1,nan_mask.shape[1]))] ,axis=0)

    indicator_mat = (padded_mat[1:,:]-padded_mat[:-1,:]).T

    pos_start = np.where(indicator_mat==1)
    pos_end = np.where(indicator_mat==-1)
    lens = (pos_end[1]-pos_start[1])[:,None]

    start_index = pos_start[1][:,None]
    time_series = pos_start[0][:,None]
    test_points = np.concatenate([start_index,time_series,lens],axis=1)
    temp = np.copy(test_points[:,2])

    if (temp.shape[0]>1):
        block_size = temp[int(temp.shape[0]/10):-int(temp.shape[0]/10)-1].mean()
    else :
        block_size = temp.mean()

    if int(block_size) <= 10:
        block_size = 10  # or some other minimum safe default
        if matrix.shape[0] <= block_size:
            block_size = 4

    w = int(10*np.log10(block_size))

    val_block_size = int(min(block_size,w))

    num_missing = int(num_missing/val_block_size)
    train_matrix = copy.deepcopy(matrix)

    val_points = []
    
    for _ in range(num_missing):
        validation_points = np.random.uniform(0, matrix.shape[0] - val_block_size, (matrix.shape[1])).astype(int)
        for i,x in enumerate(validation_points) :
            train_matrix[x:x+val_block_size,i] = np.nan
            val_points.append([x,i,val_block_size])

    return train_matrix, matrix, np.array(val_points), test_points, int(block_size), w

