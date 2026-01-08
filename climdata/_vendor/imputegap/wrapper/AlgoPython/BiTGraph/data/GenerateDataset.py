# ===============================================================================================================
# SOURCE: https://github.com/chenxiaodanhit/BiTGraph
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://openreview.net/pdf?id=O9nZCwdGcG
# ===============================================================================================================


import os
from torch.utils.data import Dataset
from torch.utils.data import DataLoader
import numpy as np
import torch

from imputegap.tools import utils

np.set_printoptions(threshold=np.inf)

class StandardScaler:
    """
    Standard the input
    """

    def __init__(self, mean, std):
        self.mean = mean
        self.std = std

    def transform(self, data):
        return (data - self.mean) / self.std

    def inverse_transform(self, data):
        return (data * self.std) + self.mean
class TSDataset(Dataset):

    def __init__(self, Data, Label=None, mask=None, masks_target=None):
        self.Data = Data
        self.Label = Label
        self.mask = mask
        self.masks_target = masks_target


    def __len__(self):
        return len(self.Data)

    def __getitem__(self, index):
        data = torch.Tensor(self.Data[index])
        if self.Label is not None:
            label = torch.Tensor(self.Label[index])
        else:
            label = None
        mask = torch.Tensor(self.mask[index])
        if self.masks_target is not None:
            masks_target = torch.Tensor(self.masks_target[index])
        else:
            masks_target = None

        return data,label,mask,masks_target


def get_0_1_array(array,rate=0.2):
    zeros_num = int(array.size * rate)
    new_array = np.ones(array.size)
    new_array[:zeros_num] = 0
    np.random.shuffle(new_array)
    re_array = new_array.reshape(array.shape)
    return re_array

def synthetic_data(mask_ratio, dataset, tr_ratio=0.9, tag="alpha", verbose=False):

    if isinstance(dataset, str):
        data_list = []
        file_path = os.path.join(os.path.dirname(__file__), 'electricity.txt')
        print(f"Validation paper: {file_path}")
        with open(file_path, 'r') as f:
            reader = f.readlines()
            for row in reader:
                data_list.append(row.split(','))

        data = np.array(data_list).astype('float')
        mask = get_0_1_array(data, mask_ratio)
        data = data[:, :, None].astype('float32')
        mask = mask[:, :, None].astype('int32')
        return data, mask
    else:
        cont_data_matrix = dataset.copy()

        nan_replacement = -999999
        artificial_training_drop = 0.30
        offset = 0.05

        # building test set ================================================================================================
        original_missing_ratio = utils.get_missing_ratio(cont_data_matrix)
        cont_data_matrix, new_mask, error = utils.prepare_testing_set(incomp_m=cont_data_matrix, original_missing_ratio=original_missing_ratio, tr_ratio=tr_ratio, verbose=False)

        valid_rows = ~np.isnan(cont_data_matrix).any(axis=1)
        data_tr = cont_data_matrix[valid_rows]
        mask_tr = new_mask[valid_rows]
        has_observed = np.any(mask_tr == 1)

        if verbose:
            print(f"{data_tr.shape = }")
            print(f"{mask_tr.shape = }")
            print("Contains observed values:", has_observed)

        gt_data_matrix = utils.prevent_leakage(cont_data_matrix, new_mask, nan_replacement, False)
        # building test set ================================================================================================

        mask_train = utils.generate_random_mask(gt=data_tr, mask_test=mask_tr, mask_valid=mask_tr, droprate=artificial_training_drop, offset=offset, verbose=False)

        data = gt_data_matrix[:, :, None].astype('float32')
        data_c = data_tr[:, :, None].astype('float32')

        new_mask = 1 - new_mask
        mask_train = 1 - mask_train

        new_mask = new_mask[:, :, None].astype('int32')
        mask_train = mask_train[:, :, None].astype('int32')

        return data, data_c, mask_train, new_mask


def split_data_by_ratio(x,y, mask,mask_target,val_ratio, test_ratio):
    idx = np.arange(x.shape[0])
    # print('idx shape:',idx.shape)
    idx_shuffle = idx.copy()
    # np.random.shuffle(idx_shuffle)
    data_len = x.shape[0]
    test_x = x[idx_shuffle[-int(data_len * test_ratio):]]
    test_y = y[idx_shuffle[-int(data_len * test_ratio):]]
    test_x_mask = mask[idx_shuffle[-int(data_len * test_ratio):]]
    test_y_mask = mask_target[idx_shuffle[-int(data_len * test_ratio):]]

    val_x = x[idx_shuffle[-int(data_len * (test_ratio + val_ratio)):-int(data_len * test_ratio)]]
    val_y = y[idx_shuffle[-int(data_len * (test_ratio + val_ratio)):-int(data_len * test_ratio)]]
    val_x_mask = mask[idx_shuffle[-int(data_len * (test_ratio + val_ratio)):-int(data_len * test_ratio)]]
    val_y_mask = mask_target[idx_shuffle[-int(data_len * (test_ratio + val_ratio)):-int(data_len * test_ratio)]]

    train_x = x[idx_shuffle[:-int(data_len * (test_ratio + val_ratio))]]
    train_y = y[idx_shuffle[:-int(data_len * (test_ratio + val_ratio))]]
    train_x_mask = mask[idx_shuffle[:-int(data_len * (test_ratio + val_ratio))]]
    train_y_mask = mask_target[idx_shuffle[:-int(data_len * (test_ratio + val_ratio))]]

    return train_x,train_y,train_x_mask,train_y_mask,val_x,val_y,val_x_mask,val_y_mask,test_x,test_y,test_x_mask,test_y_mask


def Add_Window_Horizon(data,mask, window=3, horizon=1):
    '''
    :param data: shape [B, ...]
    :param window:
    :param horizon:
    :return: X is [B, W, ...], Y is [B, H, ...]
    '''
    length = len(data)
    end_index = length - horizon - window + 1
    X = []      #windows
    Y = []  #horizon
    masks=[]
    masks_target=[]
    index = 0

    while index < end_index:
        X.append(data[index:index+window])
        masks.append(mask[index:index+window])
        Y.append(data[index+window:index+window+horizon])
        masks_target.append(mask[index+window:index+window+horizon])
        index = index + 1
    X = np.array(X)  #backcast B,W,N,D
    Y = np.array(Y)  #forecast B,H,N,D
    masks = np.array(masks)
    masks_target=np.array(masks_target)

    return X, Y, masks, masks_target


def reverse_window_horizon(X, Y, window=3, horizon=1):
    '''
    Reconstruct the original time series from overlapped windows (backcast + forecast).

    :param X: shape [B, window, ...]
    :param Y: shape [B, horizon, ...]
    :param window: size of backcast window
    :param horizon: size of forecast horizon
    :return: reconstructed_data: shape [original_length, ...]
    '''
    B = X.shape[0]
    feature_shape = X.shape[2:]
    original_length = B + window + horizon - 1

    # Initialize reconstruction arrays
    data_recon = np.zeros((original_length,) + feature_shape)
    count_recon = np.zeros((original_length,) + feature_shape)

    for i in range(B):
        # Backcast window
        data_recon[i:i + window] += X[i]
        count_recon[i:i + window] += 1

        # Forecast horizon
        data_recon[i + window:i + window + horizon] += Y[i]
        count_recon[i + window:i + window + horizon] += 1

    # Avoid division by zero
    count_recon[count_recon == 0] = 1
    reconstructed = data_recon / count_recon

    return reconstructed.squeeze()




def loaddataset_imputegap(history_len, pred_len, mask_ratio, dataset, batch_size=32, num_workers=0, tr_ratio=0.9, verbose=True, deep_verbose=False):

    if deep_verbose:
        print("\n\nLoading and transforming dataset...")

    data_imp, data_tr, mask_train_imputegap, mask_test_imputegap = synthetic_data(mask_ratio, dataset, tr_ratio=tr_ratio, tag="beta", verbose=deep_verbose)

    if deep_verbose:
        print(f"\n\t{data_imp.shape = }")
        print(f"\t{mask_test_imputegap.shape = }")
        print("\n")
        print(f"\t{data_tr.shape = }")
        print(f"\t{mask_train_imputegap.shape = }")

    num_zeros = np.sum(data_tr == 0)
    num_zeros_data = np.sum(data_imp == 0)
    if deep_verbose:
        print("\tNumber of zeros:", num_zeros)
        print("\tNumber of zeros num_zeros_data :", num_zeros_data)

    x_tr, y_tr, mask_2_tr, mask_target_2_tr = Add_Window_Horizon(data_tr, mask_train_imputegap, history_len, pred_len)
    x_ts, y_ts, mask_2_ts, mask_target_2_ts = Add_Window_Horizon(data_imp, mask_test_imputegap, history_len, pred_len)

    if verbose:
        print(f"\n\t{x_tr.shape = }")
        print(f"\t{mask_2_tr.shape = }")
        print(f"\t{x_ts.shape = }")
        print(f"\t{mask_2_ts.shape = }")

    train_x, train_y, masks_tra, masks_target_tra, val_x, val_y, masks_val, masks_target_val, test_x, test_y, masks_test, masks_target_test = split_data_by_ratio(x_tr, y_tr, mask_2_tr, mask_target_2_tr, 0, 0.3)

    if deep_verbose:
        print(f"\n\t{train_x.shape =}")
        print(f"\t{test_x.shape =}")

    if batch_size >= dataset.shape[0] // 2:
        batch_size = dataset.shape[0] // 10

    if batch_size < 1 :
        batch_size = 1

    if verbose:
        print(f"\t{batch_size = }\n")

    scaler = StandardScaler(mean=data_tr.mean(), std=data_tr.std())
    x_tra = scaler.transform(train_x)
    y_tra = scaler.transform(train_y)
    x_val = scaler.transform(val_x)
    y_val = scaler.transform(val_y)
    x_test = scaler.transform(test_x)
    y_test = scaler.transform(test_y)

    train_dataset = TSDataset(x_tra, y_tra, masks_tra, masks_target_tra)
    val_dataset = TSDataset(x_val, y_val, masks_val, masks_target_val)
    test_dataset = TSDataset(x_test, y_test, masks_test, masks_target_test)

    train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, drop_last=False)
    val_dataloader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, drop_last=False)
    test_dataloader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, num_workers=num_workers, drop_last=False)

    test_val = reverse_window_horizon(x_ts, y_ts, window=history_len, horizon=pred_len)
    if deep_verbose:
        print(f"{test_val.shape = }")
        print(f"{data_imp.shape = }")
        print(f"{np.allclose(data_imp.squeeze(), test_val, atol=1e-6) = }")


    scaler_ts = StandardScaler(mean=x_ts.mean(), std=x_ts.std())
    x_tst_imp = scaler_ts.transform(x_ts)
    y_tst_imp = scaler_ts.transform(y_ts)



    imputegap_dataset = TSDataset(x_tst_imp, y_tst_imp, mask_2_ts, mask_target_2_ts)
    imputegap_dataloader = DataLoader(imputegap_dataset, batch_size=batch_size, shuffle=False, num_workers=0, drop_last=False)

    return train_dataloader, val_dataloader, test_dataloader, imputegap_dataloader, scaler, scaler


def loaddataset(history_len,pred_len,mask_ratio,dataset,num_workers=0):
    data_numpy,mask=synthetic_data(mask_ratio,dataset)
    x, y, mask,mask_target = Add_Window_Horizon(
        data_numpy, mask,history_len, pred_len)

    print(f"{data_numpy.shape = }")
    print(f"\n\t{x.shape = }")
    print(f"\t{y.shape = }")

    train_x,train_y,masks_tra,masks_target_tra,val_x,val_y,masks_val,masks_target_val,test_x,test_y,masks_test,masks_target_test = split_data_by_ratio(x,y, mask,mask_target, 0.1, 0.2)

    print(f"{train_x.shape = }")
    print(f"\n\t{val_x.shape = }")
    print(f"\t{test_x.shape = }")

    scaler = StandardScaler(mean=train_x.mean(), std=train_x.std())
    x_tra = scaler.transform(train_x)
    y_tra = scaler.transform(train_y)
    x_val = scaler.transform(val_x)
    y_val = scaler.transform(val_y)
    x_test = scaler.transform(test_x)
    y_test = scaler.transform(test_y)

    train_dataset = TSDataset(x_tra, y_tra,masks_tra,masks_target_tra)
    val_dataset = TSDataset(x_val, y_val,masks_val,masks_target_val)
    test_dataset = TSDataset(x_test, y_test,masks_test,masks_target_test)
    train_dataloader = DataLoader(train_dataset, batch_size=32, shuffle=True, num_workers=num_workers, drop_last=True)
    val_dataloader = DataLoader(val_dataset, batch_size=32, shuffle=False, num_workers=num_workers, drop_last=True)
    test_dataloader = DataLoader(test_dataset,batch_size=32, shuffle=False, num_workers=num_workers, drop_last=False)


    return train_dataloader, val_dataloader, test_dataloader, test_dataloader, scaler, scaler


if __name__ == '__main__':
    print('')