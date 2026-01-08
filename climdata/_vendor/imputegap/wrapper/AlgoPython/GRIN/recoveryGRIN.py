# ===============================================================================================================
# SOURCE: https://github.com/Graph-Machine-Learning-Group/grin
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://openreview.net/pdf?id=kOu3-S3wJ7
# ===============================================================================================================


import datetime
import os
import pathlib
import numpy as np
import pytorch_lightning as pl
import torch
import torch.nn.functional as F
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger
from torch.optim.lr_scheduler import CosineAnnealingLR

from imputegap.recovery.manager import TimeSeries
from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.GRIN.lib.data.datamodule import SpatioTemporalDataModule
from imputegap.wrapper.AlgoPython.GRIN.lib.data.imputation_dataset import ImputationDataset, GraphImputationDataset
from imputegap.wrapper.AlgoPython.GRIN.lib.nn import models
from imputegap.wrapper.AlgoPython.GRIN.lib.nn.utils.metric_base import MaskedMetric
from imputegap.wrapper.AlgoPython.GRIN.lib.nn.utils.metrics import MaskedMAE, MaskedMAPE, MaskedMSE, MaskedMRE
from imputegap.wrapper.AlgoPython.GRIN.lib.utils import parser_utils

from imputegap.wrapper.AlgoPython.GRIN.lib import datasets, fillers, config
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def has_graph_support(model_cls):
    return model_cls in [models.GRINet]


def get_model_classes(model_str):
    if model_str == 'grin':
        model, filler = models.GRINet, fillers.GraphFiller
    else:
        raise ValueError(f'Model {model_str} not available.')
    return model, filler


def recoveryGRIN(input, d_hidden=32, lr=0.001, batch_size=-1, window=1, alpha=10.0, patience=4, epochs=20, workers=0,
                 adj_threshold=0.1, val_len=0.2, test_len=0.2, d_ff=16, ff_dropout=0.1, stride=1, l2_reg=0.0,
                 grad_clip_val=5.0, grad_clip_algorithm="norm", loss_fn="l1_loss", use_lr_schedule=True, hint_rate=0.7,
                 g_train_freq=1, d_train_freq=5, tr_ratio=0.9, seed=42, verbose=True):

    recov = np.copy(input)
    m_mask = np.isnan(input)

    if batch_size == -1:
        batch_size = utils.compute_batch_size(data=input, min_size=4, max_size=32, divisor=4, verbose=verbose)

    if verbose:
        print(f"\n(IMPUTATION) GRIN\n\tMatrix: {input.shape[0]}, {input.shape[1]}\n\tbatch_size: {batch_size}\n\tlr: {lr}\n\twindow: {window}\n\talpha: {alpha}\n\tpatience: {patience}\n\tepochs: {epochs}\n\tworkers: {workers}\n")

    nan_row_selector = np.any(np.isnan(input), axis=1)
    cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(input, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=0.0, prevent_leak=True, offset=0.05, seed=seed, verbose=False)
    if error:
        return input

    input_data = np.copy(cont_data_matrix)
    M, N = input_data.shape

    # Get indices where the value is True
    s = ~nan_row_selector
    tr_indice = np.where(s)[0]
    np.random.seed(42)

    if M <= 4:
        num_to_flip = len(tr_indice) // 2
    else:
        num_to_flip = len(tr_indice) // 4

    other_indices = np.random.choice(tr_indice, size=num_to_flip, replace=False)
    flip_indices = np.setdiff1d(tr_indice, other_indices)
    training_mask = np.zeros((M, N), dtype=np.uint8)  # or dtype=bool if preferred
    training_mask[flip_indices, :] = 1

    if window > N :
        window = 1

    torch.set_num_threads(seed)
    pl.seed_everything(seed)

    model_cls, filler_cls = get_model_classes('grin')
    dataset = datasets.MissingValuesMyData(cont_data_matrix, training_mask, mask_test)

    ########################################
    # create logdir and save configuration #
    ########################################

    # Define split configuration
    split_conf = {
        "lr": lr,
        "epochs": epochs,
        "patience": patience,
        "l2_reg": l2_reg,
        "grad_clip_val": grad_clip_val,
        "grad_clip_algorithm": grad_clip_algorithm,
        "loss_fn": loss_fn,
        "use_lr_schedule": use_lr_schedule,
        "adj_threshold": adj_threshold,
        "alpha": alpha,
        "hint_rate": hint_rate,
        "g_train_freq": g_train_freq,
        "d_train_freq": d_train_freq,
        "val_len": M,
        "test_len": len(flip_indices),
        "window": window,
        "stride": stride,
        "d_hidden": d_hidden,  # Default or replace with correct value
        "d_ff": d_ff,  # Default or replace with correct value
        "ff_dropout": ff_dropout  # Default or replace with correct value
    }

    exp_name = f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}_{seed}"
    logdir = os.path.join(config['logs'], 'grin', exp_name)
    pathlib.Path(logdir).mkdir(parents=True)

    ########################################
    # data module                          #
    ########################################

    # instantiate dataset
    dataset_cls = GraphImputationDataset if has_graph_support(model_cls) else ImputationDataset
    torch_dataset = dataset_cls(*dataset.numpy(return_idx=True), mask=dataset.training_mask, eval_mask=dataset.eval_mask, window=window, stride=stride,)

    indices = np.arange(M)
    train_idxs = flip_indices
    test_idxs = other_indices
    val_idxs = indices

    if verbose:
        print(f"🔍 torch size: {len(torch_dataset)}")
        print(f"🔍 Training Indices: {len(train_idxs) if train_idxs is not None else 0}")
        print(f"🔍 Test Indices: {len(test_idxs) if test_idxs is not None else 0}")
        print(f"🔍 Validation Indices: {len(val_idxs) if val_idxs is not None else 0}")

    if len(torch_dataset) <= 1 or len(test_idxs) <= 1 :
        print(f"\n(ERROR) Number of series to train to small for GRIN: {len(test_idxs)}\n\tPlease increase the number of series or change the dataset used.\n")
        return input


    data_module_conf = { "scale": True, "scaling_axis": "global", "scaling_type": "std", "scale_exogenous": None, "train_idxs": train_idxs, "val_idxs": val_idxs, "test_idxs": test_idxs, "batch_size": batch_size, "workers": workers, "samples_per_epoch": None, "verbose": verbose }

    dm = SpatioTemporalDataModule(torch_dataset, **data_module_conf, )
    dm.setup()

    adj = dataset.get_similarity(thr=adj_threshold)
    np.fill_diagonal(adj, 0.0)

    ########################################
    # predictor                            #
    ########################################

    # model's inputs
    additional_model_hparams = dict(adj=adj, d_in=dm.d_in, n_nodes=dm.n_nodes)
    model_kwargs = parser_utils.filter_args(args={**split_conf, **additional_model_hparams}, target_cls=model_cls  )

    # loss and metrics
    loss_fn = MaskedMetric(metric_fn=getattr(F, loss_fn), metric_kwargs={'reduction': 'none'})

    metrics = {'mae': MaskedMAE(compute_on_step=False), 'mape': MaskedMAPE(compute_on_step=False), 'mse': MaskedMSE(compute_on_step=False), 'mre': MaskedMRE(compute_on_step=False)}

    # filler's inputs
    scheduler_class = CosineAnnealingLR if use_lr_schedule else None
    additional_filler_hparams = dict(model_class=model_cls,
                                     model_kwargs=model_kwargs,
                                     optim_class=torch.optim.Adam,
                                     optim_kwargs={'lr': lr, 'weight_decay': l2_reg},
                                     loss_fn=loss_fn,
                                     metrics=metrics,
                                     scheduler_class=scheduler_class,
                                     scheduler_kwargs={'eta_min': 0.0001, 'T_max': epochs },
                                     alpha=alpha,
                                     hint_rate=hint_rate,
                                     g_train_freq=g_train_freq,
                                     d_train_freq=d_train_freq)

    filler_kwargs = parser_utils.filter_args(args={**split_conf, **additional_filler_hparams}, target_cls=filler_cls, return_dict=True)
    filler = filler_cls(**filler_kwargs)
    filler.verbose = verbose

    ########################################
    # training                             #
    ########################################
    if verbose:
        print("\ntraining...")

    ########################################

    # callbacks
    early_stop_callback = EarlyStopping(monitor='val_loss', patience=patience, mode='min')
    checkpoint_callback = ModelCheckpoint(dirpath=logdir, save_top_k=1, monitor='val_loss', mode='min')

    logger = TensorBoardLogger(logdir, name="model")

    trainer = pl.Trainer(max_epochs=epochs,
                         logger=True,
                         default_root_dir=logdir,
                         accelerator="gpu" if torch.cuda.is_available() else "cpu",  # Automatically detect GPU/CPU
                         devices=1 if torch.cuda.is_available() else "auto", # Use 1 GPU if available, otherwise default
                         gradient_clip_val=grad_clip_val,
                         gradient_clip_algorithm=grad_clip_algorithm,
                         callbacks=[early_stop_callback, checkpoint_callback],
                         enable_progress_bar=verbose,  # Disable tqdm bars
                         enable_model_summary=verbose,  # Disable model summary print
                         enable_checkpointing=True,
                         log_every_n_steps=(M%batch_size))

    trainer.fit(filler, datamodule=dm)

    ########################################
    # testing                              #
    ########################################
    if verbose:
        print("\nreconstruct...")

    checkpoint = torch.load(checkpoint_callback.best_model_path, map_location=device, weights_only=False)  # weights_only=False by default
    filler.load_state_dict(checkpoint["state_dict"])

    filler.freeze()
    trainer.test(datamodule=dm, ckpt_path="best")
    filler.eval()

    if torch.cuda.is_available():
        filler.cuda()

    with torch.no_grad():
        y_true, y_hat, mask = filler.predict_loader(dm.val_dataloader(), return_mask=True)

    # Debugging the shapes before reshaping
    y_hat = y_hat.detach().cpu().numpy().reshape(input_data.shape)

    print(f"{y_hat.shape = }")

    if verbose:
        print("🔍 y_hat shape before reshape:", y_hat.shape)
        print("🔍 Expected input_data shape:", input_data.shape)

    y_hat = y_hat

    recov[m_mask] = y_hat[m_mask]

    return recov

