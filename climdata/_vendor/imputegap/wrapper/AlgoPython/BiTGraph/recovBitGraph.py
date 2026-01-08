# ===============================================================================================================
# SOURCE: https://github.com/chenxiaodanhit/BiTGraph
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://openreview.net/pdf?id=O9nZCwdGcG
# ===============================================================================================================


import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from imputegap.wrapper.AlgoPython.BiTGraph.models.BiaTCGNet.BiaTCGNet import Model
from imputegap.wrapper.AlgoPython.BiTGraph.data.GenerateDataset import loaddataset, loaddataset_imputegap, reverse_window_horizon
from imputegap.tools import utils

torch.multiprocessing.set_sharing_strategy('file_system')
device = "cuda" if torch.cuda.is_available() else "cpu"
criterion = nn.L1Loss().to(device)


def train(model, seq_len, pred_len, mask_ratio, batch_size, epochs, imputegap_dataset=None, num_workers=0, tr_ratio=0.9, verbose=True, deep_verbose=False):

    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    torch.set_num_threads(1)

    if isinstance(imputegap_dataset, str):
        train_dataloader, val_dataloader, test_dataloader, imputegap_dataloader, scaler, scaler_ts = loaddataset(seq_len, pred_len, mask_ratio, imputegap_dataset, num_workers=num_workers)
    else:
        train_dataloader, val_dataloader, test_dataloader, imputegap_dataloader, scaler, scaler_ts = loaddataset_imputegap(seq_len, pred_len, mask_ratio, imputegap_dataset, batch_size, num_workers=num_workers, tr_ratio=tr_ratio, verbose=verbose, deep_verbose=deep_verbose)

    if verbose:
        print(f"\t{len(train_dataloader) = }")
        print(f"\t{len(test_dataloader) = }")
        print(f"\t{len(imputegap_dataloader) = }")

    best_loss, k, best_model = 9999999.99, 0, None

    if verbose:
        print("\n\ntraining...")

    for epoch in range(epochs):
        model.train()

        for i, (x, y, mask, target_mask) in enumerate(train_dataloader):
            x, y, mask, target_mask = x.to(device), y.to(device), mask.to(device), target_mask.to(device)
            x = x * mask
            y = y * target_mask
            x_hat = model(x, mask, k)
            loss = torch.sum(torch.abs(x_hat - y) * target_mask) / torch.sum(target_mask)
            optimizer.zero_grad()  # optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        loss = evaluate(model, test_dataloader, scaler)

        if verbose:
            print('\t\tepoch, loss:', epoch, loss)

        if (loss < best_loss):
            best_loss = loss
            if verbose:
                print('\t\t\tbest_loss:', best_loss)

    if isinstance(imputegap_dataset, str):
        imputed_matrix, full_y = reconstruct(model, val_dataloader, scaler_ts, "reconstruct", verbose)
    else:
        imputed_matrix, full_y = reconstruct(model, imputegap_dataloader, scaler_ts, "reconstruct", verbose)

    return imputed_matrix, full_y, model


def evaluate(model, val_iter, scaler, tag="testing"):
    model.eval()
    loss=0.0
    k=0
    with torch.no_grad():
        for i, (x, y, mask, target_mask) in enumerate(val_iter):
            x, y, mask, target_mask = x.to(device), y.to(device), mask.to(device), target_mask.to(device)
            x = x * mask
            x_hat=model(x, mask, k)
            x_hat = scaler.inverse_transform(x_hat)
            y = scaler.inverse_transform(y)
            losses = torch.sum(torch.abs(x_hat-y)*target_mask)/torch.sum(target_mask)
            loss+=losses
    return loss/len(val_iter)

def reconstruct(model, imputegap_dataloader, scaler, tag="reconstruct", verbose=True):
    model.eval()
    k = 0
    loss = 0.0
    imputed_batches = []
    y_batches = []

    if verbose:
        print(f"\n\n{tag}...")

    with torch.no_grad():
        for i, (x, y, mask, mask_y) in enumerate(imputegap_dataloader):
            x, mask = x.to(device), mask.to(device)
            y, mask_y = y.to(device), mask_y.to(device)

            x = x * mask
            y = y * mask_y

            x_hat = model(x, mask, k)
            k += 1

            x_hat = scaler.inverse_transform(x_hat)
            y = scaler.inverse_transform(y)

            losses = torch.sum(torch.abs(x_hat - y) * mask_y) / torch.sum(mask_y)
            loss += losses

            # Move to CPU and convert to numpy
            imp = x_hat.detach().cpu().numpy()
            y = y.detach().cpu().numpy()

            imputed_batches.append(imp)
            y_batches.append(y)

    full_imputed_matrix = np.concatenate(imputed_batches, axis=0)
    full_y = np.concatenate(y_batches, axis=0)

    if verbose:
        print(f'total loss: {loss:.4f}')

    return full_imputed_matrix, full_y





def recoveryBitGRAPH(input=None, node_number=-1, kernel_set=[2,3,6,7], dropout=0.05, subgraph_size=5, node_dim=3, seq_len=-1, lr=0.001, epoch=10, batch_size=32, num_workers=0, tr_ratio=0.9, seed=42, verbose=True):

    num_layers = 2
    deep_verbose = False

    test_error = input.copy()
    _, _, _, _, error = utils.dl_integration_transformation(test_error, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, seed=seed, verbose=False)
    if error:
        return input

    if isinstance(input, str):
        node_number = 321
        data = input
        pred_len = 24

        if verbose:
            print(f"(IMPUTATION) BitGRAPH\n\tMatrix:{input.shape[0]}, {input.shape[1]}\n\tnode_number: {node_number}\n\tkernel_set: {kernel_set}\n\tdropout: {dropout}\n\tsubgraph_size\n\t{subgraph_size}\n\tnode_dim: {node_dim}\n\tseq_len: {seq_len}\n\tlr: {lr}\n\tepochs: {epoch}\n\tpred_len: {pred_len}\n\tnum_workers: {num_workers}\n\rtr_ratio: {tr_ratio}\n\tseed {seed}")

    else:
        data = np.copy(input)
        recov = np.copy(input)
        missing_mask = np.isnan(input)

        if data.shape[0] <= 9:
            print(f"\n(ERROR) Number of series to train to small for LLMs: current {data.shape[0]}\n\tPlease increase the number of series or change the dataset used.\n")
            return input

        if node_number == -1:
            node_number = data.shape[1]

        enc_in = data.shape[1]
        dec_in = data.shape[1]
        c_out = data.shape[1]
        affine = 1
        pred_len = 1

        kernel_size = max(kernel_set)  # Use largest kernel for worst-case
        min_required = (kernel_size - 1) * num_layers + 1 + 1
        if min_required > input.shape[1] // 4:
            kernel_set = [1]

        if seq_len == -1:
            seq_len = 24
            if data.shape[0]//4 <= seq_len:
                kernel_size = max(kernel_set)  # Use largest kernel for worst-case
                min_required = (kernel_size - 1) * num_layers + 1 + 1
                seq_len = max(input.shape[0] // 4, min_required)

                if data.shape[0] <= 10:
                    seq_len = 2
                    kernel_set = [1]
        else:
            seq_len = seq_len


        if verbose:
            print(f"\n(IMPUTATION) BitGRAPH: Matrix Shape: ({input.shape[0]}, {input.shape[1]})"
                  f"\n\tnode_number: {node_number}\n\tkernel_set: {kernel_set}\n\tdropout: {dropout}\n\t"
                  f"subgraph_size: {subgraph_size}\n\tnode_dim: {node_dim}\n\tseq_len: {seq_len}\n\t"
                  f"lr: {lr}\n\tepochs: {epoch}\n\tpred_len: {pred_len}\n\tseed {seed}\n\t"
                  f"enc_in: {enc_in}\n\tdec_in: {dec_in}\n\tc_out: {c_out}\n\tbatch_size: {batch_size}\n\tnum_workers: {num_workers}")

    model = Model(True, True, 2, node_number, kernel_set, device=device,
                  predefined_A=None, dropout=dropout, subgraph_size=subgraph_size, node_dim=node_dim,
                  dilation_exponential=1, conv_channels=16, residual_channels=16, skip_channels=16, end_channels=32,
                  seq_length=seq_len, in_dim=1, out_len=pred_len, out_dim=1, layers=num_layers, propalpha=0.05,
                  tanhalpha=3, layer_norm_affline=affine)  # 2 4 6

    model.to(device)

    imputed_matrix, full_y, best_model = train(model, seq_len, pred_len, 0.2, batch_size, epoch, data, num_workers=num_workers, tr_ratio=tr_ratio, verbose=verbose, deep_verbose=deep_verbose)

    if verbose:
        print("Imputed Matrix Shape Before Reshaping:", imputed_matrix.shape)

    reshaped_imputed_matrix = reverse_window_horizon(imputed_matrix, full_y, window=seq_len, horizon=pred_len)

    if verbose:
        print("Reshaped Imputed Matrix Shape:", reshaped_imputed_matrix.shape)

    recov[missing_mask] = reshaped_imputed_matrix[missing_mask]
    assert np.allclose(recov[~missing_mask], input[~missing_mask]), "Observed values were changed!"

    return recov


