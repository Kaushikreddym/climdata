# ===============================================================================================================
# SOURCE: https://github.com/jsyoon0823/MRNN
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://ieeexplore.ieee.org/document/8485748
# ===============================================================================================================

import numpy as np
from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.MRNN import Data_Loader
from imputegap.wrapper.AlgoPython.MRNN.M_RNN import M_RNN

def reconstruct(recover, x, m, size, seq_length, verbose=True):

    if verbose:
        print(f"\n{recover.shape = }")
        print(f"{x.shape = }")
        print(f"{m = }")
        print(f"{size = }")
        print(f"{seq_length = }")

    # part 1: upper block
    for si in range(0, seq_length - 1):  # si = sequence index
        i = size + si  # index in the main matrix
        for j in range(0, m):
            if np.isnan(x[i][j]):
                val = 0.0
                for sj in range(0, si + 1):
                    val += recover[sj][si - sj][j]
                x[i][j] = val / (si + 1)  # average

    # part 2: middle block
    for ri in range(seq_length - 1, len(recover)):  # ):
        i = size + ri

        if i >= x.shape[0]:
            i = x.shape[0]-1

        for j in range(0, m):
            if np.isnan(x[i][j]):
                val = 0.0
                for sj in range(0, seq_length):
                    val += recover[ri - sj][sj][j]
                x[i][j] = val / seq_length  # average

    # part 3: lower block
    for si in range(0, seq_length):  # si = sequence index
        i = len(x) - si - 1  # index in the main matrix
        ri = len(recover) - si - 1
        for j in range(0, m):
            if np.isnan(x[i][j]):
                val = 0.0
                for sj in range(0, si + 1):
                    val += recover[ri + sj][seq_length - sj - 1][j]
                x[i][j] = val / (si + 1)  # average

    return np.asarray(x)

def mrnn_recov(matrix_in, hidden_dim=10, learning_rate=0.01, iterations=1000, seq_length=7, tr_ratio=0.9, verbose=True, seed=42):

    recov = np.copy(matrix_in)
    m_mask = np.isnan(matrix_in)

    x = False
    deep_verbose = verbose
    while not x:
        # ==================================================================================================================
        cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(matrix_in, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, seed=seed, verbose=False)
        # ==================================================================================================================
        if error:
            return matrix_in

        nan_row_selector = np.any(np.isnan(cont_data_matrix), axis=1)
        cont_data_test = cont_data_matrix[nan_row_selector]
        cont_data_train = cont_data_matrix[~nan_row_selector]

        if cont_data_train.shape[0] <= 2:
            print(f"\n(ERROR) Number of series to train to small for MRNN: {cont_data_train.shape[0]}\n\tPlease increase the number of series or change the dataset used.\n")
            return matrix_in


        indices_ts = np.where(nan_row_selector)[0]

        if hidden_dim == -1:
            hidden_dim = cont_data_train.shape[1]
        if seq_length > mask_train.shape[0]//3:
            seq_length = 7

            if seq_length > cont_data_train.shape[0]//3:
                seq_length = cont_data_train.shape[0]//4
            if seq_length <= 1:
                seq_length = 2

        if deep_verbose:
            print(f"(IMPUTATION) MRNN\n\tMatrix: {matrix_in.shape[0]}, {matrix_in.shape[1]}\n\thidden_dim: {hidden_dim}\n\tlearning_rate: {learning_rate}\n\titerations: {iterations}\n\tseq_length: {seq_length}")
            print(f"\n{cont_data_train.shape = }")
            print(f"{cont_data_test.shape = }\n")

        #_, trainZ, trainM, trainT, testX, testZ, testM, testT, train_size, test_size
        trainX, trainZ, trainM, trainT, testX, testZ, testM, testT, train_size, test_size = Data_Loader.Data_Loader_With_Dataset(seq_length, cont_data_train, cont_data_test, mask_train, nan_val=None, verbose=verbose)


        if len(testX) > 0:
            x = True
        else:
            tr_ratio = tr_ratio - 0.1
            print(f"Test set is empty. As training ratio is to high, we reduce its value to tr_ratio={tr_ratio}\n")
            deep_verbose = False
        if tr_ratio == 0:
            break

    Recover_trainX, Recover_testX = M_RNN(trainZ, trainM, trainT, testZ, testM, testT, hidden_dim=hidden_dim, learning_rate=learning_rate, iterations=iterations)

    m = len(recov[1])
    alpha = recov.copy()
    beta = recov.copy()

    if verbose:
        print("\nreconstruction...")

    perfect_size = min(train_size, test_size)
    tr = reconstruct(Recover_trainX, alpha, m, 0, seq_length, verbose=False)
    ts = reconstruct(Recover_testX, alpha, m, perfect_size, seq_length, verbose=False)

    tr_inc, ts_inc = 0, 0
    for idx, row in enumerate(recov):
        if idx in indices_ts:
            beta[idx] = ts[ts_inc]
            ts_inc += 1
        else:
            beta[idx] = tr[tr_inc]
            tr_inc += 1

    recov[m_mask] = beta[m_mask]

    return recov
