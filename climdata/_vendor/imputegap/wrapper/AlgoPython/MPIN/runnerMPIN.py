# ===============================================================================================================
# SOURCE: https://github.com/XLI-2020/MPIN
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://www.vldb.org/pvldb/vol17/p345-li.pdf
# ===============================================================================================================

import pandas as pd
import torch

from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.MPIN.utils.load_dataset import get_model_size
from imputegap.wrapper.AlgoPython.MPIN.utils.regressor import MLPNet
from imputegap.wrapper.AlgoPython.MPIN.utils.DynamicGNN import DynamicGCN, DynamicGAT, DynamicGraphSAGE, StaticGCN, StaticGraphSAGE, StaticGAT

import torch.optim as optim
import copy
import numpy as np
import random


def knn_graph_torch(x, k, loop=False):
    device = x.device  # ensure all tensors are on the same device

    num_nodes = x.size(0)
    dist = torch.cdist(x, x, p=2)  # shape: [N, N]

    if not loop:
        dist.fill_diagonal_(float('inf'))

    knn_idx = dist.topk(k, largest=False).indices  # shape: [N, k]

    row = torch.arange(num_nodes, device=device).unsqueeze(1).repeat(1, k).flatten()
    col = knn_idx.flatten()
    edge_index = torch.stack([row, col], dim=0)  # shape: [2, N*k]

    return edge_index


def data_transform(X, X_mask, eval_ratio=0.05):
    eval_mask = np.zeros(X_mask.shape)
    rows, cols = np.where(X_mask == 1)
    eval_row_index_index = random.sample(range(len(rows)), int(eval_ratio * len(rows)))
    eval_row_index = rows[eval_row_index_index]
    eval_col_index = cols[eval_row_index_index]
    X_mask[eval_row_index, eval_col_index] = 0
    eval_mask[eval_row_index, eval_col_index] = 1
    eval_X = copy.copy(X)
    X[eval_row_index, eval_col_index] = 0
    return X, X_mask, eval_X, eval_mask


def build_GNN(in_channels, out_channels, k, base, device):
    if base == 'GAT':
        gnn = DynamicGAT(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    elif base == 'GCN':
        gnn = DynamicGCN(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    elif base == 'SAGE':
        gnn = DynamicGraphSAGE(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    return gnn


def build_GNN_static(in_channels, out_channels, k, base, device):
    if base == 'GAT':
        gnn = StaticGAT(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    elif base == 'GCN':
        gnn = StaticGCN(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    elif base == 'SAGE':
        gnn = StaticGraphSAGE(in_channels=in_channels, out_channels=out_channels, k=k).to(device)
    return gnn


def get_window_data(base_X, base_X_mask, start, end, ratio):
    X = base_X[int(len(base_X) * start * ratio):int(len(base_X) * end * ratio)]
    X_mask = base_X_mask[int(len(base_X) * start * ratio):int(len(base_X) * end * ratio)]
    return X, X_mask


def window_imputation(input, input_test, mask, mask_eval, start, end, sample_ratio, initial_state_dict=None, X_last=None, mask_last=None,
                      mae_last=None, transfer=False, lr=0.01, weight_decay=0.1, epochs=200, out_channels=256,
                      state=True, k=10, base="SAGE", thre=0.25, eval_ratio=0.05, dynamic=False, device=None, verbose=True):

    X, X_mask = get_window_data(base_X=input, base_X_mask=mask, start=start, end=end, ratio=sample_ratio)

    ori_X = copy.copy(X)
    feature_dim = ori_X.shape[1]
    ori_X_mask = copy.copy(X_mask)

    all_mask = copy.copy(X_mask)
    all_X = copy.copy(X)

    if X_last:
        X_last = np.array(X_last)
        all_X = np.concatenate([X_last, X], axis=0)
        all_mask = np.concatenate([mask_last, X_mask], axis=0)

        X_last = X_last.tolist()

    all_mask_ts = torch.FloatTensor(all_mask).to(device)

    gram_matrix = torch.mm(all_mask_ts, all_mask_ts.t())  # compute the gram product

    gram_vec = gram_matrix.diagonal()

    gram_row_sum = gram_matrix.sum(dim=0)

    value_vec = gram_vec - (gram_row_sum - gram_vec) / (gram_matrix.shape[0] - 1)

    keep_index = torch.where(value_vec > thre * (feature_dim - 1))[0]
    keep_index = keep_index.data.cpu().numpy()

    keep_mask = all_mask[keep_index]

    keep_X = all_X[keep_index]

    X, X_mask, eval_X, eval_mask = data_transform(X, X_mask, eval_ratio=eval_ratio)
    #eval_mask = get_window_data(base_X=input, base_X_mask=mask_eval, start=start, end=end, ratio=sample_ratio)


    if X_last:
        X_last = np.array(X_last)
        shp_last = X_last.shape
        eval_X = np.concatenate([X_last, eval_X], axis=0)
        X = np.concatenate([X_last, X], axis=0)
        eval_mask_last = np.zeros(shp_last)
        eval_mask = np.concatenate([eval_mask_last, eval_mask], axis=0)
        X_mask = np.concatenate([mask_last, X_mask], axis=0)

    in_channels = X.shape[1]
    out_channels = X.shape[1]
    X = torch.FloatTensor(X).to(device)
    X_mask = torch.LongTensor(X_mask).to(device)
    eval_X = torch.FloatTensor(eval_X).to(device)
    eval_mask = torch.LongTensor(eval_mask).to(device)

    # build model
    if dynamic == 'true':
        gnn = build_GNN(in_channels=in_channels, out_channels=out_channels, k=k, base=base, device=device)
        gnn2 = build_GNN(in_channels=in_channels, out_channels=out_channels, k=k, base=base, device=device)
    else:
        gnn = build_GNN_static(in_channels=in_channels, out_channels=out_channels, k=k, base=base, device=device)
        gnn2 = build_GNN_static(in_channels=in_channels, out_channels=out_channels, k=k, base=base, device=device)

    model_list = [gnn, gnn2]
    regressor = MLPNet(out_channels, in_channels).to(device)

    if initial_state_dict != None:
        gnn.load_state_dict(initial_state_dict['gnn'])
        gnn2.load_state_dict(initial_state_dict['gnn2'])
        if not transfer:
            regressor.load_state_dict(initial_state_dict['regressor'])

    trainable_parameters = []
    for model in model_list:
        trainable_parameters.extend(list(model.parameters()))

    trainable_parameters.extend(list(regressor.parameters()))
    filter_fn = list(filter(lambda p: p.requires_grad, trainable_parameters))

    num_of_params = sum(p.numel() for p in filter_fn)

    model_size = get_model_size(gnn) + get_model_size(gnn2) + get_model_size(regressor)
    model_size = round(model_size, 6)

    num_of_params = num_of_params / 1e6

    opt = optim.Adam(filter_fn, lr=lr, weight_decay=weight_decay)

    graph_impute_layers = len(model_list)

    X_knn = copy.deepcopy(X)

    edge_index = knn_graph_torch(X_knn, k, loop=False)

    min_mae_error = 1e9
    min_mse_error = None
    min_mape_error = None
    opt_epoch = None
    opt_time = None
    best_X_imputed = None
    best_state_dict = None

    for pre_epoch in range(epochs):
        gnn.train()
        gnn2.train()
        regressor.train()
        opt.zero_grad()
        loss = 0
        X_imputed = copy.copy(X)

        for i in range(graph_impute_layers):
            if dynamic == 'true':
                X_emb = model_list[i](X_imputed)
            else:
                X_emb, edge_index = model_list[i](X_imputed, edge_index)

            pred = regressor(X_emb)
            X_imputed = X * X_mask + pred * (1 - X_mask)
            temp_loss = torch.sum(torch.abs(X - pred) * X_mask) / (torch.sum(X_mask) + 1e-5)
            loss += temp_loss

        loss.backward()
        opt.step()
        train_loss = loss.item()
        if verbose:
            print('\n\t\t\t{n} epoch loss:'.format(n=pre_epoch), train_loss, '.............')

        trans_X = copy.copy(X_imputed)
        trans_eval_X = copy.copy(eval_X)

        epoch_state_dict = {'gnn': gnn.state_dict(), 'gnn2': gnn2.state_dict(), 'regressor': regressor.state_dict()}

        gnn.eval()
        gnn2.eval()
        regressor.eval()

        with torch.no_grad():

            mae_error = torch.sum(torch.abs(trans_X - trans_eval_X) * eval_mask) / torch.sum(eval_mask)
            mse_error = torch.sum(((trans_X - trans_eval_X) ** 2) * eval_mask) / torch.sum(eval_mask)
            mape_error = torch.sum(torch.abs(trans_X - trans_eval_X) * eval_mask ) / (torch.sum(torch.abs(trans_eval_X) * eval_mask) + 1e-12)
            if verbose:
                print('\t\t\t\timputegap impute error MAE:', mae_error.item())
                print('\t\t\t\timputegap impute error MSE:', mse_error.item())
                print('\t\t\t\timputegap impute error MRE:', mape_error.item())


            if mae_error.item() < min_mae_error:
                opt_epoch = copy.copy(pre_epoch)
                min_mae_error = round(mae_error.item(), 6)
                if verbose:
                    print('\t\t\t{epoch}_opt_mae_error'.format(epoch=pre_epoch), min_mae_error)

                min_mse_error = round(mse_error.item(), 6)
                min_mape_error = round(mape_error.item(), 6)

                if verbose:
                    print('\t\t\t{epoch}_opt time:'.format(epoch=pre_epoch), opt_time)

                best_X_imputed = copy.copy(X_imputed)
                best_X_imputed = best_X_imputed.detach().cpu().numpy()
                best_X_imputed = best_X_imputed * (1 - ori_X_mask) + ori_X * ori_X_mask
                best_state_dict = copy.copy(epoch_state_dict)

    results_list = [opt_epoch, min_mae_error, min_mse_error, min_mape_error, num_of_params, model_size, opt_time, 0]

    if mae_last and (min_mae_error > mae_last) and (state == 'true'):
        best_state_dict = copy.copy(initial_state_dict)
    return best_state_dict, keep_X.tolist(), keep_mask, results_list, min_mae_error, best_X_imputed


def recoverMPIN(input, mode="alone", window=2, k=10, lr=0.01, weight_decay=0.1, epochs=200, num_of_iteration=5, thre=0.25,
                base="SAGE", out_channels=64, eval_ratio=0.05, state=True, dynamic=True, tr_ratio=0.9, seed=0, verbose=True):

    recov = np.copy(input)
    m_mask = np.isnan(input)
    full_imputed = np.zeros_like(input)
    full_mask = np.zeros_like(input, dtype=bool)  # to track what has been filled

    if verbose:
        print(f"(IMPUTATION) MPIN\n\tMatrix: {input.shape[0]}, {input.shape[1]}\n\tk: {k}\n\twindow: {window}\n\tlr: {lr}\n\tweight: {weight_decay}\n\tbase: {base}\n\tthreshold: {thre}\n\tepochs: {epochs}")

    torch.random.manual_seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    random.seed(seed)
    cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(input, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=-1, prevent_leak=True, offset=0.05, seed=seed, verbose=False)
    if error:
        return input

    num_windows = window

    results_schema = ['opt_epoch', 'opt_mae', 'mse', 'mape', 'para', 'memo', 'opt_time', 'tot_time']

    iter_results_list = []

    for iteration in range(num_of_iteration):
        results_collect = []
        for w in range(num_windows):
            if w == 0:
                window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test, start=w, end=w + 1, sample_ratio=1 / num_windows, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)
            else:
                if mode == 'alone':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test, start=w, end=w + 1, sample_ratio=1 / num_windows, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

                elif mode == 'data':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test,start=w, end=w + 1, sample_ratio=1 / num_windows, X_last=X_last, mask_last=mask_last, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

                elif mode == 'state':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test,start=w, end=w + 1, sample_ratio=1 / num_windows, initial_state_dict=window_best_state, mae_last=mae_last, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

                elif mode == 'state+transfer':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test,start=w, end=w + 1, sample_ratio=1 / num_windows, initial_state_dict=window_best_state, transfer=True, mae_last=mae_last, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

                elif mode == 'data+state':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test,start=w, end=w + 1, sample_ratio=1 / num_windows, initial_state_dict=window_best_state, X_last=X_last, mask_last=mask_last, mae_last=mae_last, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels, state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

                elif mode == 'data+state+transfer':
                    window_best_state, X_last, mask_last, window_results, mae_last, best_X = window_imputation(input=cont_data_matrix, input_test=cont_data_matrix, mask=mask_train, mask_eval=mask_test,start=w, end=w + 1, sample_ratio=1 / num_windows, initial_state_dict=window_best_state, X_last=X_last, mask_last=mask_last, transfer=True, mae_last=mae_last, lr=lr, weight_decay=weight_decay, epochs=epochs, out_channels=out_channels,state=state, k=k, base=base, thre=thre, eval_ratio=eval_ratio, dynamic=dynamic, device=device, verbose=verbose)

            results_collect.append(window_results)

            if window > 1:
                start_idx = int(len(input) * w / num_windows)
                end_idx = int(len(input) * (w + 1) / num_windows)
                full_imputed[start_idx:end_idx] = best_X
                full_mask[start_idx:end_idx] = True
            else:
                full_imputed = best_X

        df = pd.DataFrame(results_collect, index=range(num_windows), columns=results_schema)
        iter_results_list.append(df)

    if verbose:
     print('\nreconstruct...')

    full_imputed = np.array(full_imputed)
    print(f"{full_imputed.shape = }")

    recov[m_mask] = full_imputed[m_mask]

    return recov
