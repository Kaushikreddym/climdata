# ===============================================================================================================
# SOURCE: https://github.com/xuangu-fang/BayOTIDE
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/abs/2308.14906
# ===============================================================================================================


import os
import numpy as np
import torch 
import sys
import imputegap
from imputegap.recovery.manager import TimeSeries
from imputegap.tools import utils
sys.path.append("../")
import tqdm
import yaml
torch.random.manual_seed(300)
from imputegap.wrapper.AlgoPython.BayOTIDE import utils_BayOTIDE, model_BayOTIDE
import warnings
warnings.filterwarnings("ignore")
torch.random.manual_seed(300)


def generate_mask(data_matrix, drop_rate=0.8, valid_rate=0.2, verbose=False):
    """
    Generate train/test/valid masks based only on existing NaN positions in data_matrix.
    """
    nan_mask = np.isnan(data_matrix)
    nan_indices = np.argwhere(nan_mask)
    np.random.shuffle(nan_indices)

    n = len(nan_indices)
    n_test = int(n * drop_rate)
    n_valid = int(n * valid_rate)
    n_train = n - n_test - n_valid

    if verbose:
        print(f"\n{n =}")

    # Sanity check
    if verbose:
        print(f"\n{n_test =}")
        print(f"{n_valid =}")
        print(f"{n_train =}\n")

    train_idx = nan_indices[:n_train]
    test_idx = nan_indices[n_train:n_train + n_test]
    valid_idx = nan_indices[n_train + n_test:]

    if verbose:
        print(f"\n{train_idx.shape =}")
        print(f"{test_idx.shape =}")
        print(f"{valid_idx.shape =}\n")

    mask_train = np.zeros_like(data_matrix, dtype=np.uint8)
    mask_test = np.zeros_like(data_matrix, dtype=np.uint8)
    mask_valid = np.zeros_like(data_matrix, dtype=np.uint8)

    mask_train[tuple(train_idx.T)] = 1
    mask_test[tuple(test_idx.T)] = 1
    mask_valid[tuple(valid_idx.T)] = 1

    # Sanity check
    if verbose:
        print(f"\nTrain mask NaNs: {mask_train.sum()}")
        print(f"Test mask NaNs: {mask_test.sum()}")
        print(f"Valid mask NaNs: {mask_valid.sum()}")

        print(f"{mask_train.shape =}")
        print(f"{mask_test.shape =}")
        print(f"{mask_valid.shape =}\n")

    return mask_train, mask_test, mask_valid


def recovBayOTIDE(incomp_m, K_trend=None, K_season=None, n_season=None, K_bias=None, time_scale=None, a0=None, b0=None, v=None, num_workers=0, tr_ratio=0.9, config=None, args=None, verbose=True):
    """
    Run BayOTIDE model using a provided NumPy data matrix instead of loading from a file.

    :param data_matrix: Preloaded NumPy matrix containing time series data (N x T).
    :param K_trend: Number of trend factors (optional, overrides config if provided).
    :param K_season: Number of seasonal factors (optional, overrides config if provided).
    :param n_season: Number of seasonal components per factor (optional, overrides config if provided).
    :param K_bias: Number of bias factors (optional, overrides config if provided).
    :param time_scale: Scaling factor for the time step (optional, overrides config if provided).
    :param a0: Prior hyperparameter for variance scaling (optional, overrides config if provided).
    :param b0: Prior hyperparameter for variance scaling (optional, overrides config if provided).
    :param v: Variance hyperparameter for noise modeling (optional, overrides config if provided).
    :param tr_ratio: ratio of the training set.
    :param config: Dictionary containing hyperparameters (optional).
    :param args: Parsed arguments for the model (optional).
    :return: Imputed time series matrix (N x T).
    """

    final_result = incomp_m.copy()
    mask_original_nan = np.isnan(incomp_m)

    nan_replacement             = -99999
    artificial_training_drop    = 0.4
    ts_ratio                    = 0.9
    val_ratio                   = 1-ts_ratio
    offset                      = 0.05

    if verbose:
        print(f"(IMPUTATION) BayOTIDE\n\tMatrix: {incomp_m.shape[0]}, {incomp_m.shape[1]}\n\tK_trend: {K_trend}\n\tK_season: {K_season}\n\tn_season: {n_season}\n\tK_bias: {K_bias}\n\ttime_scale: {time_scale}\n\ta0: {a0}\n\tb0: {b0}\n\tnum_workers: {num_workers}\n\ttr_ratio: {tr_ratio}\n")

    # building test set ================================================================================================
    gt_data_matrix = incomp_m.copy()
    cont_data_matrix = incomp_m.copy()

    original_missing_ratio = utils.get_missing_ratio(cont_data_matrix)
    cont_data_matrix, new_mask, error = utils.prepare_testing_set(incomp_m=cont_data_matrix, original_missing_ratio=original_missing_ratio, tr_ratio=tr_ratio, verbose=False)
    if error:
        return incomp_m
    gt_data_matrix = utils.prevent_leakage(gt_data_matrix, new_mask, nan_replacement, False)
    # building test set ================================================================================================

    sub_tensor = torch.from_numpy(gt_data_matrix).float()
    zero_ratio = (sub_tensor == nan_replacement).sum().item() / sub_tensor.numel()

    if config is None:
        # Get directory of current file
        base_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(base_dir, "config_bayotide.yaml")

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        if verbose:
            print(f"\nloading of the configuration file : {config_path}")

    # Conditional updates
    if K_trend is not None:
        config["K_trend"] = K_trend
    if K_season is not None:
        config["K_season"] = K_season
    if n_season is not None:
        config["n_season"] = n_season
    if K_bias is not None:
        config["K_bias"] = K_bias
    if time_scale is not None:
        config["time_scale"] = time_scale
    if a0 is not None:
        config["a0"] = a0
    if b0 is not None:
        config["b0"] = b0
    if v is not None:
        config["v"] = v

    device = "cuda" if torch.cuda.is_available() else "cpu"
    if verbose:
        print(f"Using device: {device}\n")

    # Set device in config
    config["device"] = device
    seed = config["seed"]

    data_save = {}
    data_save['ndims'] = sub_tensor.shape
    data_save['raw_data'] = sub_tensor
    data_save['data'] = []
    data_save['time_uni'] = np.linspace(0, 1, sub_tensor.shape[1])

    for i in range(config["num_fold"]):
        # building masks================================================================================================
        mask_test, mask_valid, nbr_nans = utils.split_mask_bwt_test_valid(cont_data_matrix, test_rate=ts_ratio, valid_rate=val_ratio, verbose=False, seed=seed)
        mask_train = utils.generate_random_mask(gt=gt_data_matrix, mask_test=mask_test, mask_valid=mask_valid, droprate=artificial_training_drop, offset=offset, verbose=False, seed=seed)
        # building masks================================================================================================

        data_save['data'].append({'mask_train': mask_train, 'mask_test': mask_test, 'mask_valid': mask_valid})


    base_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(base_dir, config["data_path"])
    np.save(file_path, data_save)
    config["data_path"] = file_path

    if verbose:
        print(f"Data saved to {file_path}")

    data_file = file_path
    hyper_dict = config

    INNER_ITER = hyper_dict["INNER_ITER"]
    EVALU_T = hyper_dict["EVALU_T"]

    if verbose:
        print("\n\ntraining...")

    for fold_id in range(config["num_fold"]):

        data_dict = utils_BayOTIDE.make_data_dict(hyper_dict, data_file, fold=fold_id)

        model = model_BayOTIDE.BayTIDE(hyper_dict, data_dict)

        model.reset()

        # one-pass along the time axis
        for T_id in tqdm.tqdm(range(model.T)):
            model.filter_predict(T_id)
            model.msg_llk_init()

            if model.mask_train[:, T_id].sum() > 0:  # at least one obseved data at current step
                for inner_it in range(INNER_ITER):
                    flag = (inner_it == (INNER_ITER - 1))

                    model.msg_approx_U(T_id)
                    model.filter_update(T_id, flag)

                    model.msg_approx_W(T_id)
                    model.post_update_W(T_id)

                model.msg_approx_tau(T_id)
                model.post_update_tau(T_id)

            else:
                model.filter_update_fake(T_id)

            if T_id % EVALU_T == 0 or T_id == model.T - 1:
                _, loss_dict = model.model_test(T_id)

                if verbose:
                    print("T_id = {}, train_rmse = {:.3f}, test_rmse= {:.3f}".format(T_id, loss_dict["train_RMSE"], loss_dict["test_RMSE"]))

    if verbose:
        print('\nsmoothing back...')

    model.smooth()
    model.post_update_U_after_smooth(0)

    # Run model test and get predictions
    pred, loss_dict = model.model_test(T_id)

    if verbose:
        print(f"\t{pred.shape =}")

    # Fill NaNs in original data
    if isinstance(pred, torch.Tensor):
        pred = pred.detach().cpu().numpy()

    final_result[mask_original_nan] = pred[mask_original_nan]

    if verbose:
        print(f"\t{final_result.shape =}")

    return final_result


if __name__ == "__main__":

    ts = TimeSeries()
    ts.load_series(imputegap.tools.utils.search_path("chlorine"))
    print(f"{ts.data.shape = }")

    # contaminate the time series with MCAR pattern
    ts_m = ts.Contamination.mcar(ts.data)

    imputed_data = recovBayOTIDE(ts.data, ts_m)

    from imputegap.recovery.imputation import Imputation

    imputer = Imputation.DeepLearning.BayOTIDE(ts_m)
    imputer.recov_data = imputed_data
    imputer.incomp_data = ts_m

    # compute and print the imputation metrics
    imputer.score(ts.data, imputed_data)
    ts.print_results(imputer.metrics)

    # plot the recovered time series
    ts.plot(input_data=ts.data, incomp_data=ts_m, recov_data=imputed_data, nbr_series=9, subplot=True, algorithm=imputer.algorithm, save_path="./imputegap_assets/imputation")
    ts.plot(input_data=ts.data, nbr_series=9,  save_path="./imputegap_assets/imputation")
    ts.plot(input_data=imputed_data, nbr_series=9,  save_path="./imputegap_assets/imputation")