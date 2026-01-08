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

def mean_recovery(matrix):
    means = np.nanmean(matrix,axis=0)
    for i in range(matrix.shape[1]):
        matrix[:,i] = np.nan_to_num(matrix[:,i],nan=means[i])
    return matrix


def zero_recovery(matrix):
    return np.nan_to_num(matrix)
