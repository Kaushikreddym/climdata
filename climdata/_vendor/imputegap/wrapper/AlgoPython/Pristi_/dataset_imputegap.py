# ===============================================================================================================
# SOURCE: https://github.com/LMZZML/PriSTI
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://ieeexplore.ieee.org/document/10184808
# ===============================================================================================================

from torch.utils.data import DataLoader, Dataset
import numpy as np
import torch
import torchcde

from imputegap.tools import utils

import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
device = "cuda" if torch.cuda.is_available() else "cpu"


def sample_mask(shape, p=0.0015, p_noise=0.05, max_seq=1, min_seq=1, rng=None):
    if rng is None:
        rand = np.random.random
        randint = np.random.randint
    else:
        rand = rng.random
        randint = rng.integers
    mask = rand(shape) < p
    for col in range(mask.shape[1]):
        idxs = np.flatnonzero(mask[:, col])
        if not len(idxs):
            continue
        fault_len = min_seq
        if max_seq > min_seq:
            fault_len = fault_len + int(randint(max_seq - min_seq))
        idxs_ext = np.concatenate([np.arange(i, i + fault_len) for i in idxs])
        idxs = np.unique(idxs_ext)
        idxs = np.clip(idxs, 0, shape[0] - 1)
        mask[idxs, col] = True
    mask = mask | (rand(mask.shape) < p_noise)
    return mask.astype('uint8')


class ImputeGAP_Dataset(Dataset):
    def __init__(self, data, tr_ratio=0.9, eval_length=1, mode="train", val_len=0.1, test_len=0.2, missing_pattern='block', is_interpolate=False, target_strategy='random', batch_size=1):

        cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(data, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=0.9, split_val=0.1, nan_val=-999999, prevent_leak=True, offset=0.05, seed=42, verbose=False)


        self.eval_length = eval_length
        self.is_interpolate = is_interpolate
        self.target_strategy = target_strategy
        self.mode = mode
        self.target = cont_data_matrix.shape[1]

        train_values = cont_data_matrix[mask_train.astype(bool)]

        if train_values.size > 0:
            self.train_mean = np.array([np.nanmean(train_values)])  # Wrap in 1D array
            self.train_std = np.array([np.nanstd(train_values)])  # Wrap in 1D array
        else:
            self.train_mean = np.array([0.0])
            self.train_std = np.array([1.0])  # Avoid division by zero

        # create data for batch
        self.use_index = []
        self.cut_length = []
        self.data = data
        self.batch_size = batch_size

        SEED = 9101112

        if mode == 'train':
            self.observed_mask = ~mask_train
            self.gt_mask = mask_train
            self.observed_data = cont_data_matrix
        elif mode == 'valid':
            self.observed_mask = ~mask_train
            self.gt_mask = mask_test
            self.observed_data = cont_data_matrix
        elif mode == 'test':
            self.observed_mask = ~mask_train
            self.gt_mask = mask_val
            self.observed_data = cont_data_matrix
        else:
            self.observed_mask = ~mask_train
            self.gt_mask = mask_test+mask_val
            self.observed_data = cont_data_matrix

        current_length = len(self.observed_mask) - eval_length + 1

        if mode == "test":
            n_sample = len(self.observed_data) // eval_length
            c_index = np.arange(0, 0 + eval_length * n_sample, eval_length)
            self.use_index += c_index.tolist()
            self.cut_length += [0] * len(c_index)
            if len(self.observed_data) % eval_length != 0:
                self.use_index += [current_length - 1]
                self.cut_length += [eval_length - len(self.observed_data) % eval_length]
        elif mode != "test":
            self.use_index = np.arange(current_length)
            self.cut_length = [0] * len(self.use_index)

    def __getitem__(self, org_index):
        index = self.use_index[org_index]
        ob_data = self.observed_data[index: index + self.eval_length]
        ob_mask = self.observed_mask[index: index + self.eval_length]
        gt_mask = self.gt_mask[index: index + self.eval_length]
        cond_mask = torch.tensor(gt_mask).to(torch.float32)
        timepoints= np.arange(self.eval_length)
        cut_length = self.cut_length[org_index]

        #ob_data = self.observed_data
        #ob_mask = self.observed_mask
        #gt_mask = self.gt_mask
        #cond_mask = torch.tensor(self.gt_mask).to(torch.float32)
        #timepoints= np.arange(self.batch_size)
        #cut_length = self.cut_length[0]*self.data.shape[0]

        s = {
            "observed_data": ob_data,
            "observed_mask": ob_mask,
            "gt_mask": gt_mask,
            "timepoints": timepoints,
            "cut_length": cut_length,
            "cond_mask": cond_mask,
            "index": self.use_index
        }

        if self.is_interpolate:
            tmp_data = torch.tensor(ob_data).to(torch.float64)
            itp_data = torch.where(cond_mask == 0, float('nan'), tmp_data).to(torch.float32)
            itp_data = torchcde.linear_interpolation_coeffs(itp_data.permute(1, 0).unsqueeze(-1)).squeeze(-1).permute(1, 0)
            s["coeffs"] = itp_data.numpy()
        return s

    def __len__(self):
        return len(self.use_index)


def get_dataloader(data, tr_ratio, batch_size, device, val_len=0.1, test_len=0.2, missing_pattern='block', is_interpolate=False, num_workers=0, target_strategy='random', eval_length=1):
    dataset = ImputeGAP_Dataset(data, tr_ratio=tr_ratio, eval_length=eval_length, mode="train", val_len=val_len, test_len=test_len, missing_pattern=missing_pattern, is_interpolate=is_interpolate, target_strategy=target_strategy, batch_size=batch_size)
    train_loader = DataLoader(dataset, batch_size=batch_size, num_workers=num_workers, shuffle=False)

    dataset_test = ImputeGAP_Dataset(data, tr_ratio=tr_ratio, eval_length=eval_length,  mode="test", val_len=val_len, test_len=test_len, missing_pattern=missing_pattern, is_interpolate=is_interpolate, target_strategy=target_strategy, batch_size=batch_size)
    test_loader = DataLoader(dataset_test, batch_size=batch_size, num_workers=num_workers, shuffle=False)

    dataset_valid = ImputeGAP_Dataset(data, tr_ratio=tr_ratio, eval_length=eval_length, mode="valid", val_len=val_len, test_len=test_len, missing_pattern=missing_pattern, is_interpolate=is_interpolate, target_strategy=target_strategy, batch_size=batch_size)
    valid_loader = DataLoader(dataset_valid, batch_size=batch_size, num_workers=num_workers, shuffle=False)

    dataset_recon = ImputeGAP_Dataset(data, tr_ratio=tr_ratio, eval_length=eval_length, mode="reconstruct", val_len=val_len, test_len=test_len, missing_pattern=missing_pattern, is_interpolate=is_interpolate, target_strategy=target_strategy, batch_size=batch_size)
    recon_loader = DataLoader(dataset_recon, batch_size=batch_size, num_workers=num_workers, shuffle=False)

    scaler = torch.from_numpy(dataset.train_std).to(device).float()

    mean_scaler = torch.from_numpy(dataset.train_mean).to(device).float()

    return train_loader, valid_loader, test_loader, recon_loader, scaler, mean_scaler, dataset.target

