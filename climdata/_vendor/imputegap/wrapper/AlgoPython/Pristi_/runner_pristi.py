# ===============================================================================================================
# SOURCE: https://github.com/LMZZML/PriSTI
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://ieeexplore.ieee.org/document/10184808
# ===============================================================================================================

import logging
import torch
import datetime
import json
import yaml
import os
import numpy as np

from imputegap.wrapper.AlgoPython.Pristi_.dataset_imputegap import get_dataloader
from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.Pristi_.main_model import PriSTI_PemsBAY
from imputegap.wrapper.AlgoPython.Pristi_.utils import train, reconstruct


def recov_pristi(data, target_strategy="block", unconditional=True, batch_size=-1, embedding=-1, num_workers=0, nsample=100, modelfolder="", tr_ratio=0.9, seed=42, verbose=True):

    recov = np.copy(data)
    m_mask = np.isnan(data)

    test_error = data.copy()
    _, _, _, _, error = utils.dl_integration_transformation(test_error, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, seed=seed, verbose=False)
    if error:
        return data

    if batch_size == -1:
        batch_size = utils.compute_batch_size(data, min_size=2, max_size=16, verbose=verbose)
    if embedding == -1:
        embedding = utils.compute_batch_size(data, min_size=4, max_size=32, verbose=False)
        if embedding %2 == 1:
            embedding = embedding + 1

    if data.shape[0] % batch_size != 0:
        batch_size = utils.compute_batch_size(data, min_size=2, max_size=10, verbose=False)
        if data.shape[0] % batch_size != 0:
            batch_size = 1

    eval_length = data.shape[0] // batch_size

    if verbose:
        print(f"(IMPUTATION) PriSTI\n\tMatrix Shape: {data.shape[0]}, {data.shape[1]}\n\ttarget_strategy: {target_strategy}\n\tunconditional: {unconditional}\n\tnum_workers: {num_workers}\n\tbatch_size: {batch_size}\n\tembedding: {embedding}\n\teval_length: {eval_length}\n\ttr_ratio: {tr_ratio}\n\tseed: {seed}\n")

    SEED = seed
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    torch.cuda.manual_seed(SEED)

    base_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_dir, "config", "traffic.yaml")
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    config["model"]["is_unconditional"] = unconditional
    config["model"]["target_strategy"] = target_strategy
    config["diffusion"]["adj_file"] = 'pems-bay'
    config["seed"] = SEED
    device = "cuda" if torch.cuda.is_available() else "cpu"
    config["train"]["batch_size"] = batch_size
    config["diffusion"]["diffusion_embedding_dim"] = embedding
    config["model"]["timeemb"] = embedding
    config["model"]["featureemb"] = embedding//2

    #print(json.dumps(config, indent=4))

    current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    foldername = ("./imputegap_assets/models/pristi_" + current_time + "/")
    os.makedirs(foldername, exist_ok=True)
    with open(foldername + "config.json", "w") as f:
        json.dump(config, f, indent=4)

    train_loader, valid_loader, test_loader, recon_loader, scaler, mean_scaler, target = get_dataloader(data, tr_ratio=tr_ratio, batch_size=config["train"]["batch_size"], device=device, missing_pattern="block", is_interpolate=config["model"]["use_guide"], num_workers=num_workers, target_strategy=target_strategy, eval_length=eval_length)

    model = PriSTI_PemsBAY(train_loader=train_loader, config=config, device=device, target_dim=target).to(device)

    if modelfolder == "":
        train(model, config["train"], train_loader, valid_loader=valid_loader, foldername=foldername, )
    else:
        model.load_state_dict(torch.load("./imputegap_assets/models/" + modelfolder + "/model.pth"))

    logging.basicConfig(filename=foldername + '/test_model.log', level=logging.DEBUG)
    logging.info("model_name={}".format(modelfolder))

    if verbose:
        print("\n\nreconstruction")

    reconstruction = reconstruct(model, data, recon_loader, batch_size, nsample=1, device=device)

    if verbose:
        print(f"\n{reconstruction.shape = }")

    recov[m_mask] = reconstruction[m_mask]

    return recov