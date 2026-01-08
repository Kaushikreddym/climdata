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

import time
import warnings
warnings.filterwarnings("ignore")
args = utils_BayOTIDE.parse_args_dynamic_streaming()

torch.random.manual_seed(args.seed)

from imputegap.recovery.imputation import Imputation


def validation_random_mask(shape, drop_rate=0.2, valid_rate=0.1):
    """
    train_ratio: 1-valid_rate-drop_rate
    test_ratio: drop_rate
    valid_ratio: valid_rate
    """
    N, T = shape

    np.random.seed(300)

    mask_train_list = []
    mask_test_list = []
    mask_valid_list = []

    for t in range(T):
        mask = np.random.rand(N)
        mask_train = np.where(mask > drop_rate + valid_rate, 1, 0)
        mask_test = np.where(mask < drop_rate, 1, 0)
        mask_valid = np.where((mask > drop_rate) & (mask < drop_rate + valid_rate), 1, 0)

        mask_train_list.append(mask_train)
        mask_test_list.append(mask_test)
        mask_valid_list.append(mask_valid)

    mask_train = np.stack(mask_train_list, axis=1)
    mask_test = np.stack(mask_test_list, axis=1)
    mask_valid = np.stack(mask_valid_list, axis=1)

    return mask_train, mask_test, mask_valid



def valBayOTIDE(ts, gt, data_path, K_trend=None, K_season=None, n_season=None, K_bias=None, time_scale=None, a0=None, b0=None, v=None, tr_ratio=0.6, config=None, args=None, verbose=True):
    """
    validation of BayOTIDE model using a provided NumPy data matrix instead of loading from a file.

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



    if data_path == "guangzhou_impute_r_0.2.npy" or data_path == "simu_impute_r_0.1.npy":
        base_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(base_dir, data_path)
        full_data = np.load(data_path, allow_pickle=True).item()
        gt = full_data["raw_data"]
        ts.data = gt.copy()

    gt_data_matrix = gt.copy()
    cont_data_matrix= gt.copy()
    final_result = gt.copy()

    sub_tensor = torch.from_numpy(gt_data_matrix).float()

    if config is None:
        # Get directory of current file
        base_dir = os.path.dirname(os.path.abspath(__file__))
        config_path = os.path.join(base_dir, "config_bayotide.yaml")

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        if verbose:
            print(f"\n3.a) loading of the configuration file : {config_path = }")

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
        print(f"3.b) Using device: {device}\n")

    # Set device in config
    config["device"] = device
    seed = config["seed"]

    data_save = {}
    data_save['ndims'] = sub_tensor.shape

    data_save['data'] = []
    data_save['time_uni'] = np.linspace(0, 1, sub_tensor.shape[1])

    for i in range(config["num_fold"]):
        mask_train, mask_test, mask_valid = validation_random_mask(sub_tensor.shape, 0.4, 0.1)
        data_save['data'].append({'mask_train': mask_train, 'mask_test': mask_test, 'mask_valid': mask_valid})

        print("\nNumber of 1s in mask_train:", np.sum(mask_train))
        print("Number of 1s in mask_test:", np.sum(mask_test))
        print("Number of 1s in mask_valid:", np.sum(mask_valid), "\n")

    sub_tensor[mask_test == 1] = 0
    sub_tensor[mask_valid == 1] = 0

    data_save['raw_data'] = sub_tensor

    if data_path is None:
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
                    print("T_id = {}, train_rmse = {:.3f}, test_rmse= {:.3f}".format(T_id, loss_dict["train_RMSE"],
                                                                                     loss_dict["test_RMSE"]))

    if verbose:
        print('\n\n5.a) smoothing back...')

    model.smooth()
    model.post_update_U_after_smooth(0)

    # Run model test and get predictions
    pred, loss_dict = model.model_test(T_id)

    if verbose:
        print(f"{pred.shape =}")

    if isinstance(pred, torch.Tensor):
        pred = pred.detach().cpu().numpy()

    # Fill NaNs in original data
    final_result[mask_test] = pred[mask_test]
    cont_data_matrix[mask_test == 1] = np.nan

    if verbose:
        print(f"{final_result.shape =}")

    imputer = Imputation.DeepLearning.BayOTIDE(cont_data_matrix)
    imputer.recov_data = final_result

    imputer.score(gt, final_result)

    ts.print_results(imputer.metrics)
    ts.plot(ts.data, cont_data_matrix, final_result, subplot=True, nbr_series=9)

    print("end______________")

    return final_result


if __name__ == "__main__":

    ts = TimeSeries()

    imputed_data = valBayOTIDE(ts, gt=None, data_path="guangzhou_impute_r_0.2.npy")
    imputed_data = valBayOTIDE(ts, gt=None, data_path="simu_impute_r_0.1.npy")

    ts.load_series(imputegap.tools.utils.search_path("chlorine"))
    print(f"{ts.data.shape = }")
    imputed_data = valBayOTIDE(ts, ts.data, data_path=None)


