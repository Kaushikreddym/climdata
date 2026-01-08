# ===============================================================================================================
# SOURCE: https://github.com/Chengyui/NuwaTS/tree/master
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/pdf/2405.15317
# ===============================================================================================================

from imputegap.wrapper.AlgoPython.NuwaTS.dataprovider.data_loader_imputation import Dataset_Custom
from torch.utils.data import DataLoader

data_dict = {'custom': Dataset_Custom}

def data_provider(args, flag, tr=None, ts=None, m_tr=None, m_ts=None, ts_m=None, verbose=False):
    Data = data_dict[args.data]
    timeenc = 0 if args.embed != 'timeF' else 1

    shuffle_flag = False
    drop_last = False
    batch_size = args.batch_size
    freq = args.freq

    data_set = Data(
        root_path=args.root_path,
        data_path=args.data_path,
        flag=flag,
        size=[args.seq_len, args.label_len, args.pred_len, args.patch_size],
        features=args.features,
        target=args.target,
        timeenc=timeenc,
        freq=freq,
        ts_m=ts_m,
        batch_size=batch_size,
        tr=tr, ts=ts, m_tr=m_tr, m_ts=m_ts, verbose=verbose
    )

    data_loader = DataLoader(
            data_set,
            batch_size=batch_size,
            shuffle=shuffle_flag,
            num_workers=args.num_workers,
            drop_last=drop_last)

    return data_set, data_loader
