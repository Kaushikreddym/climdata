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

from imputegap.recovery.imputation import Imputation
from imputegap.wrapper.AlgoPython.DeepMVI.transformer import transformer_recovery
from imputegap.recovery.manager import TimeSeries
from imputegap.tools import utils


def deep_mvi_recovery(input, max_epoch=1000, patience=2, lr=1e-3, tr_ratio=0.9, seed=42, verbose=True):

    if verbose:
        print(f"(IMPUTATION) DEEP-MVI\n\tMatrix: {input.shape[0]}, {input.shape[1]}\n\tmax_epoch: {max_epoch}\n\tpatience: {patience}\n\tlr: {lr}\n\ttr_ratio: {tr_ratio}\n\tseed: {seed}")

    recov = np.copy(input)
    m_mask = np.isnan(input)

    cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(input, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, seed=seed, verbose=False)
    if error:
        return input

    temp = transformer_recovery(cont_data_matrix, max_epoch, patience, lr, verbose=verbose)

    recov[m_mask] = temp[m_mask]

    return recov

# end function

if __name__ == '__main__' :
    ts = TimeSeries()
    ts.load_series(utils.search_path("chlorine"))
    print(f"{ts.data.shape = }")

    # contaminate the time series with MCAR pattern
    ts_m = ts.Contamination.mcar(ts.data)

    imputer = Imputation.DeepLearning.DeepMVI(ts_m)
    imputation_m = deep_mvi_recovery(ts_m)
    imputer.recov_data = imputation_m

    # compute and print the imputation metrics
    imputer.score(ts.data, imputer.recov_data)
    ts.print_results(imputer.metrics)

    # plot the recovered time series
    ts.plot(input_data=ts.data, incomp_data=ts_m, recov_data=imputer.recov_data, nbr_series=12, subplot=True,
            algorithm=imputer.algorithm, save_path="./imputegap_assets/imputation")

