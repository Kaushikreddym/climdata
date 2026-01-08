# ===============================================================================================================
# SOURCE: https://github.com/caow13/BRITS
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://papers.nips.cc/paper_files/paper/2018/hash/734e6bfcd358e25ac1db0a4241b95651-Abstract.html
# ===============================================================================================================


import torch
import torch.optim as optim

import numpy as np
import imputegap.wrapper.AlgoPython.BRITS.utils as utilsX
import imputegap.wrapper.AlgoPython.BRITS.models as models
import imputegap.wrapper.AlgoPython.BRITS.data_loader as data_loader
from imputegap.tools import utils

from imputegap.wrapper.AlgoPython.BRITS.data_prep_tf import prepare_dat


def train(model, input, batch_size, epochs, num_workers=0, verbose=True):
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    data_iter = data_loader.get_loader(input, batch_size=batch_size, num_workers=num_workers)
    for epoch in range(0, epochs):
        model.train()
        run_loss = 0.0

        for idx, data in enumerate(data_iter):  # 4
            data = utilsX.to_var(data)
            ret = model.run_on_batch(data, optimizer)
            run_loss += ret['loss'].data

            forward = data["forward"]
            values = forward['values']

            if verbose:
                print('\r Progress epoch {}, {:.2f}%, batch {} [{}], average loss {}'.format(epoch, (idx + 1) * 100.0 / len(data_iter), idx, values.shape, run_loss / (idx + 1.0)))

    return (model, data_iter)


def evaluate(model, val_iter):
    model.eval()
    imputations = []

    for idx, data in enumerate(val_iter):
        data = utilsX.to_var(data)
        ret = model.run_on_batch(data, None)
        imputation = ret['imputations'].data.cpu().numpy()
        imputations += imputation.tolist()

    imputations = np.asarray(imputations)
    return imputations


def brits_recovery(incomp_data, model="brits_i_univ", epoch=10, batch_size=7, nbr_features=1, hidden_layers=64, seq_length=36, num_workers=0, tr_ratio=0.9, seed=42, verbose=True):
    recov = np.copy(incomp_data)
    m_mask = np.isnan(incomp_data)

    if batch_size == -1:
        batch_size = utils.compute_batch_size(data=incomp_data, min_size=2, max_size=24, divisor=4, verbose=verbose)

    # ==================================================================================================================
    cont_data_matrix, mask_train, mask_test, mask_valid, error = utils.dl_integration_transformation(incomp_data, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.2, split_ts=1, split_val=0, nan_val=-99999, prevent_leak=-99999, offset=0.05, seed=seed, verbose=False)
    # ==================================================================================================================
    if error:
        return incomp_data

    prepare_dat(cont_data_matrix, "./imputegap_assets/models/brits.tmp", mask_train, mask_test, mask_valid)

    if incomp_data.ndim == 2 and nbr_features != 1:
        print(f"\n(ERROR) The number of features set is not correct for the dimension of the data {incomp_data.ndim} must be higher then 2\n\tNumber of feature set to 1.\n")
        nbr_features = 1


    if verbose:
        print(f"(IMPUTATION) {model.upper()}\n\tMatrix: {incomp_data.shape[0]}, {incomp_data.shape[1]}\n\tepoch: {epoch}\n\tbatch_size: {batch_size}\n\tnbr_features: {nbr_features}\n\tseq_length: {seq_length}\n\thidden_layers: {hidden_layers}\n\tnum_workers: {num_workers}\n\ttr_ratio: {tr_ratio}\n")

    model = getattr(models, model).Model(batch_size, nbr_features, hidden_layers, seq_length)

    if torch.cuda.is_available():
        model = model.cuda()

    (model, data_iter) = train(model, "./imputegap_assets/models/brits.tmp", batch_size, epoch, num_workers, verbose)

    res = evaluate(model, data_iter)

    recovery = np.squeeze(np.array(res))

    if verbose:
        print("recovery", recovery.shape)

    recov[m_mask] = recovery[m_mask]

    return recov
