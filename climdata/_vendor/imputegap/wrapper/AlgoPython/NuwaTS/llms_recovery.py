# ===============================================================================================================
# SOURCE: https://github.com/Chengyui/NuwaTS/tree/master
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/pdf/2405.15317
# ===============================================================================================================

import argparse
import sys
from imputegap.tools import utils
from matplotlib import pyplot as plt
from imputegap.wrapper.AlgoPython.NuwaTS.exp.exp_imputation import Exp_Imputation
import random
import numpy as np
import torch.multiprocessing

torch.multiprocessing.set_sharing_strategy('file_system')


def llms_recov(ts_m, seq_length=-1, patch_size=-1, batch_size=-1, pred_length=-1, label_length=-1, enc_in=10, dec_in=10, c_out=10, gpt_layers=6, num_workers=0, tr_ratio=0.9, model="NuwaTS", seed=42, verbose=True):
    recov = np.copy(ts_m)
    m_mask = np.isnan(ts_m)
    miss = np.copy(ts_m)

    cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(miss, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.5, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, block_selection=True, seed=seed, verbose=False)
    if error:
        return ts_m

    mask_train = 1 - mask_train
    mask_test = 1 - mask_test

    nan_row_selector = np.any(np.isnan(cont_data_matrix), axis=1)
    cont_data_test = cont_data_matrix[nan_row_selector]
    cont_mask_test = mask_test[nan_row_selector]
    cont_data_train = cont_data_matrix[~nan_row_selector]
    cont_mask_train = mask_train[~nan_row_selector]

    M, N = cont_data_train.shape
    if M <= 2:
        print(f"\n(ERROR) Number of series to train to small for LLMs: {M}\n\tPlease increase the number of series or change the dataset used.\n")
        return ts_m

    if seq_length == -1:
        seq_length = utils.compute_seq_length(M)
    if batch_size == -1:
        batch_size = utils.compute_batch_size(ts_m, 4, 16, 2, verbose)
    if patch_size == -1:
        for p in reversed(range(2, seq_length -1)):
            if seq_length % p == 0:
                patch_size = p
                break
        else:
            patch_size = 1
        if model != "NuwaTS":
            patch_size = 1

    if pred_length == -1:
        pred_length = (M//2) - seq_length + 1 - (M//seq_length)
        if pred_length < 1:
            pred_length = 1
    if label_length == -1:
        if seq_length > pred_length:
            label_length = seq_length - pred_length
        else:
            label_length = pred_length - seq_length
        if label_length < 1:
            label_length = 1
    if c_out == -1:
        if model == "NuwaTS":
            c_out = miss.shape[1]
        else:
            c_out = miss.shape[1] // patch_size
    if enc_in == -1:
        if model == "NuwaTS":
            enc_in = miss.shape[1]
        else:
            enc_in = miss.shape[1] // patch_size
    if dec_in == -1:
        if model == "NuwaTS":
            dec_in = miss.shape[1]
        else:
            dec_in = miss.shape[1] // patch_size



    sys.argv += [
        '--task_name', 'imputation',
        '--is_training', '1',
        '--root_path', 'imputegap',
        '--data_path', 'imputegap',
        '--model', str(model),
        '--data', 'custom',
        '--features', 'M',
        '--seq_len', str(seq_length),
        '--label_len', str(label_length),
        '--pred_len', str(pred_length),
        '--enc_in', str(enc_in),
        '--dec_in', str(dec_in),
        '--c_out', str(c_out),
        '--num_workers', str(num_workers),
        '--gpt_layers', str(gpt_layers),
        '--batch_size', str(batch_size),
        '--d_model', '768',
        '--patch_size', str(patch_size),
        '--des', 'NuwaTS_ECL',
        '--mlp', '1',
        '--learning_rate', '0.001',
        '--prefix_length', '1',
        '--checkpoints', './imputegap_assets/models/checkpoints/'
        '--prefix_tuning',
        '--cov_prompt',
    ]

    fix_seed = seed
    random.seed(fix_seed)
    torch.manual_seed(fix_seed)
    np.random.seed(fix_seed)

    parser = argparse.ArgumentParser(description=model)

    # basic config
    parser.add_argument('--task_name', type=str, default='denoise', help='task name, options:[long_term_forecast, short_term_forecast, imputation, classification, anomaly_detection, denoise]')
    parser.add_argument('--is_training', type=int, default=1, help='status')
    parser.add_argument('--model_id', type=str, default='test', help='model id')
    parser.add_argument('--model', type=str, required=True, default=model, help='model name, options: [Autoformer, Transformer, TimesNet]')

    # data loader
    parser.add_argument('--data', type=str, required=True, default='custom', help='dataset type')
    parser.add_argument('--root_path', type=str, default='./dataset/', help='root path of the data file')
    parser.add_argument('--data_path', type=str, default='electricity.csv', help='data file')
    parser.add_argument('--features', type=str, default='M', help='forecasting task, options:[M, S, MS]; M:multivariate predict multivariate, S:univariate predict univariate, MS:multivariate predict univariate')
    parser.add_argument('--target', type=str, default='OT', help='target feature in S or MS task')
    parser.add_argument('--freq', type=str, default='h', help='freq for time features encoding, options:[s:secondly, t:minutely, h:hourly, d:daily, b:business days, w:weekly, m:monthly], you can also use more detailed freq like 15min or 3h')
    parser.add_argument('--checkpoints', type=str, default='./imputegap_assets/models/checkpoints/', help='location of model checkpoints')

    # forecasting task
    parser.add_argument('--seq_len', type=int, default=96, help='input sequence length')
    parser.add_argument('--label_len', type=int, default=0, help='start token length')
    parser.add_argument('--pred_len', type=int, default=0, help='prediction sequence length')
    parser.add_argument('--seasonal_patterns', type=str, default='Monthly', help='subset for M4')

    # imputation task
    parser.add_argument('--test_mask_rate', type=float, default=0.8, help='test mask ratio')
    parser.add_argument('--max_iterations', type=int, default=10, help='max iterations')
    parser.add_argument('--max_optimization_iterations', type=int, default=10, help='max optimization iterations')
    parser.add_argument('--regularization_weight', type=float, default=0.05, help='regularization weight')

    # anomaly detection task
    parser.add_argument('--anomaly_ratio', type=float, default=0.25, help='prior anomaly ratio (%)')

    # model define
    parser.add_argument('--top_k', type=int, default=5, help='for TimesBlock')
    parser.add_argument('--num_kernels', type=int, default=6, help='for Inception')
    parser.add_argument('--enc_in', type=int, default=107, help='encoder input size')
    parser.add_argument('--dec_in', type=int, default=107, help='decoder input size')
    parser.add_argument('--c_out', type=int, default=7, help='output size')
    parser.add_argument('--d_model', type=int, default=512, help='dimension of model')
    parser.add_argument('--n_heads', type=int, default=8, help='num of heads')
    parser.add_argument('--e_layers', type=int, default=2, help='num of encoder layers')
    parser.add_argument('--d_layers', type=int, default=1, help='num of decoder layers')
    parser.add_argument('--d_ff', type=int, default=2048, help='dimension of fcn')
    parser.add_argument('--moving_avg', type=int, default=25, help='window size of moving average')
    parser.add_argument('--factor', type=int, default=1, help='attn factor')
    parser.add_argument('--distil', action='store_false', help='whether to use distilling in encoder, using this argument means not using distilling', default=True)
    parser.add_argument('--dropout', type=float, default=0.1, help='dropout')
    parser.add_argument('--embed', type=str, default='timeF', help='time features encoding, options:[timeF, fixed, learned]')
    parser.add_argument('--activation', type=str, default='gelu', help='activation')
    parser.add_argument('--output_attention', action='store_true', help='whether to output attention in ecoder')

    # optimization
    parser.add_argument('--num_workers', type=int, default=0, help='data loader num workers')
    parser.add_argument('--itr', type=int, default=1, help='experiments times')
    parser.add_argument('--train_epochs', type=int, default=10, help='train epochs')
    parser.add_argument('--batch_size', type=int, default=32, help='batch size of train input data')
    parser.add_argument('--patience', type=int, default=3, help='early stopping patience')
    parser.add_argument('--learning_rate', type=float, default=0.001, help='optimizer learning rate')
    parser.add_argument('--des', type=str, default='test', help='exp description')
    parser.add_argument('--loss', type=str, default='MSE', help='loss function')
    parser.add_argument('--lradj', type=str, default='type1', help='adjust learning rate')
    parser.add_argument('--use_amp', action='store_true', help='use automatic mixed precision training', default=False)

    # GPU
    parser.add_argument('--use_gpu', type=bool, default=True, help='use gpu')
    parser.add_argument('--gpu', type=int, default=0, help='gpu')
    parser.add_argument('--use_multi_gpu', action='store_true', help='use multiple gpus', default=False)
    parser.add_argument('--devices', type=str, default='0,1,2,3', help='device ids of multile gpus')

    # de-stationary projector params
    parser.add_argument('--p_hidden_dims', type=int, nargs='+', default=[128, 128], help='hidden layer dimensions of projector (List)')
    parser.add_argument('--p_hidden_layers', type=int, default=2, help='number of hidden layers in projector')
    # patching
    parser.add_argument('--patch_size', type=int, default=1)
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--gpt_layers', type=int, default=6)
    parser.add_argument('--ln', type=int, default=0)
    parser.add_argument('--mlp', type=int, default=0)
    parser.add_argument('--weight', type=float, default=0)
    parser.add_argument('--percent', type=int, default=5)

    #prefix tuning
    parser.add_argument('--prefix_tuning', action='store_true', help='', default=False)
    parser.add_argument('--prefix_tuningv2', action='store_true', help='', default=False)
    parser.add_argument('--continue_tuning', action='store_true', help='', default=False)
    parser.add_argument('--continue_tuningv2', action='store_true', help='', default=False)

    parser.add_argument('--frozen_lm', action='store_true', help='', default=False)
    parser.add_argument('--prefix_length',type=int, default=1)
    parser.add_argument('--train_all_lm', action='store_true', help='', default=False)
    parser.add_argument('--use_llama', action='store_true', help='', default=False)
    parser.add_argument('--use_bert', action='store_true', help='', default=False)
    parser.add_argument('--alignment', action='store_true', help='', default=False)

    #contrastive
    parser.add_argument('--con_weight', type=float, default=0.01, help='')
    parser.add_argument('--patch_con', action='store_true', help='', default=False)
    parser.add_argument('--temporal_con', action='store_true', help='', default=False)
    parser.add_argument('--flatten_con', action='store_true', help='', default=False)
    parser.add_argument('--best_con_num', type=int, default=128)
    # output learnable token
    parser.add_argument('--seq_token',type=int, default=0)
    parser.add_argument('--word_prompt', action='store_true', help='', default=False)
    parser.add_argument('--cov_prompt', action='store_true', help='', default=False)
    parser.add_argument('--output_token', action='store_true', help='', default=False)

    # test
    parser.add_argument('--test_all',action='store_true', help='', default=False)

    #forecasting
    parser.add_argument('--is_forecasting', action='store_true', help='', default=False)
    parser.add_argument('--auto_regressive', action='store_true', help='', default=False)

    parser.add_argument('--origin_missrate', type=float, default=0, help='')

    args, _ = parser.parse_known_args()
    args.use_gpu = True if torch.cuda.is_available() and args.use_gpu else False

    if args.use_gpu and args.use_multi_gpu:
        args.devices = args.devices.replace(' ', '')
        device_ids = args.devices.split(',')
        args.device_ids = [int(id_) for id_ in device_ids]
        args.gpu = args.device_ids[0]


    if verbose:
        print(f"(IMPUTATION) {model} (LLMs)\n\tMatrix: {miss.shape[0]}, {miss.shape[1]}\n\tseq_length: {seq_length}\n\tpatch_size: {patch_size}\n\tbatch_size: {batch_size}\n\tpred_length: {pred_length}\n\tlabel_length: {label_length}\n\tenc_in: {enc_in}\n\tdec_in: {dec_in}\n\tc_out: {c_out}\n\tgpt_layers: {gpt_layers}\n\tnum_workers: {num_workers}\n\ttr_ratio: {tr_ratio}\n\tseed: {seed}\n\tverbose: {verbose}\n\tGPU: {args.use_gpu}")

    Exp= Exp_Imputation

    if args.is_training:

        if verbose:
            print(f"\ntraining of the LLMs...\n")

        for ii in range(1):
            setting = '{}_{}_{}'.format(args.model, args.data, args.des, ii)

            exp = Exp(args)  # set experiments
            exp.train(setting, tr=cont_data_train, ts=None, m_tr=cont_mask_train, m_ts=None, model_name=model, verbose=verbose)

            if verbose:
                print(f"\n\nreconstruction...\n")
            pred, _, _  = exp.test(setting, tr=None, ts=cont_data_matrix, m_tr=None, m_ts=mask_test, model_name=model, verbose=verbose)
            torch.cuda.empty_cache()
    """
    else:
        if verbose:
            
        if args.test_all:
            ii = 0
            data_path = []
            data_type = []
            des = []
            for i in range(0,10):
                args.data_path = data_path[i]
                args.data = data_type[i]
                ddes = args.des+des[i]
                setting = '{}_{}_{}'.format(args.model, args.data, ddes, ii)

                exp = Exp(args)  # set experiments
                print('>>>>>>>testing : {}<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(setting))
                pred, _, _  = exp.test(setting, test=1, tr=cont_data_matrix, ts=cont_data_matrix, m_tr=mask_train, m_ts=mask_test, verbose=verbose)
            torch.cuda.empty_cache()
        else:
            ii = 0

            setting = '{}_{}_{}'.format(args.model, args.data, args.des, ii)

            exp = Exp(args)  # set experiments
            print('>>>>>>>testing : {}<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'.format(setting))
            pred, _, _ = exp.test(setting, test=1, tr=cont_data_matrix, ts=cont_data_matrix, m_tr=mask_train, m_ts=mask_test, verbose=verbose)
    """
    plt.close('all')

    # Check for NaNs
    if verbose:
        if np.isnan(pred).any():
            print("[WARNING] pred contains NaNs!")
            nan_count = np.isnan(pred).sum()
            print(f"Total NaNs in pred: {nan_count}")
            nan_locs = np.argwhere(np.isnan(pred))
            print(f"First few NaN locations (up to 10):\n{nan_locs[:10]}")
        else:
            print("[INFO] pred does not contain any NaNs.")

    #imputation = pred.reshape(-1, pred.shape[-1])

    total_sensors = miss.shape[0]
    total_timesteps = miss.shape[1]

    stride = 1
    window_size = seq_length

    imputation = np.zeros((total_sensors, total_timesteps))
    count = np.zeros((total_sensors, total_timesteps))

    for i in range(pred.shape[0]):
        start = i * stride
        end = start + window_size

        if end > total_sensors:
            # Trim the window if it exceeds the total length
            valid_len = total_sensors - start
            imputation[start:total_sensors] += pred[i][:valid_len]
            count[start:total_sensors] += 1
        else:
            imputation[start:end] += pred[i]
            count[start:end] += 1

    # Avoid division by zero
    count[count == 0] = 1
    imputation_llms = imputation / count

    if verbose:
        print(f"{pred.shape = }")
        print(f"{imputation_llms.shape = }")

    recov[m_mask] = imputation_llms[m_mask]

    return recov


