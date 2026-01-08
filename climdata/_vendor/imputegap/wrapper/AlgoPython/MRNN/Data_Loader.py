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

def Data_Loader_With_Dataset(seq_length, data_tr, data_ts, mask_train=None, nan_val=None, verbose=False):
    def MinMaxScaler(data):
        dmin = np.nanmin(data, 0)
        dmax = np.nanmax(data, 0)
        numerator = data - dmin
        denominator = dmax - dmin
        return numerator / (denominator + 1e-8), dmin, dmax

    def make_sequences(data, nan_val=None):
        col_no = data.shape[1]
        dataX, dataZ, dataM, dataT = [], [], [], []
        row_no = len(data) - seq_length

        for i in range(row_no):
            _x = data[i:i + seq_length]
            dataX.append(_x)

            # Mask (1 if observed, 0 if NaN)
            m = np.ones([seq_length, col_no])
            if nan_val is None:
                m[np.isnan(_x)] = 0
            else:
                m[_x == nan_val] = 0
            dataM.append(m)

            # Zero-masked version
            z = np.copy(_x)
            if nan_val is None:
                z[np.isnan(_x)] = 0
            else:
                z[_x == nan_val] = 0
            dataZ.append(z)

            # Time gap
            t = np.ones([seq_length, col_no])
            for j in range(col_no):
                for k in range(1, seq_length):
                    if m[k, j] == 0:
                        t[k, j] = t[k - 1, j] + 1
            dataT.append(t)

        return (
            np.array(dataX),
            np.array(dataZ),
            np.array(dataM),
            np.array(dataT)
        )

    if mask_train is not None:
        data_tr[mask_train] = np.nan

    # Process both training and testing datasets
    trainX, trainZ, trainM, trainT = make_sequences(data_tr, nan_val=nan_val)
    testX, testZ, testM, testT = make_sequences(data_ts, nan_val=nan_val)

    if verbose:
        print(f"{trainX.shape = }")
        print(f"{testX.shape = }\n")

    return [trainX, trainZ, trainM, trainT, testX, testZ, testM, testT, len(trainX), len(testX)]
