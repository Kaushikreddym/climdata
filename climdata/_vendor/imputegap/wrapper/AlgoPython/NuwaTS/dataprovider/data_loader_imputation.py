# ===============================================================================================================
# SOURCE: https://github.com/Chengyui/NuwaTS/tree/master
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/pdf/2405.15317
# ===============================================================================================================

from torch.utils.data import Dataset
from imputegap.wrapper.AlgoPython.NuwaTS.utils.timefeatures import time_features
import warnings

warnings.filterwarnings('ignore')

class Dataset_Custom(Dataset):
    def __init__(self, root_path, flag='train', size=None, features='S', data_path='ETTh1.csv', target='OT', scale=True, timeenc=0, freq='h',  seasonal_patterns=None, percent=10,  train_sensors=None, val_sensors=None, test_sensors=None, tr=None, ts=None, m_tr=None, m_ts=None, ts_m=None, batch_size=None, verbose=True):
        # size [seq_len, label_len, pred_len]
        # info

        self.seq_len = size[0]
        self.label_len = size[1]
        self.pred_len = size[2]
        self.patch_size = size[3]

        # init
        assert flag in ['train', 'test', 'val']
        type_map = {'train': 0, 'val': 1, 'test': 2}
        self.set_type = type_map[flag]

        self.features = features
        self.scale = scale
        self.timeenc = timeenc
        self.freq = freq

        self.percent = percent
        self.root_path = root_path
        self.data_path = data_path
        self.ts_m = ts_m
        self.tr = tr
        self.ts = ts
        self.m_tr = m_tr
        self.m_ts = m_ts
        self.m_tr_heal = None
        self.m_tr_cont = None
        self.verbose = verbose

        self.__read_data__()

    def __read_data__(self):
        import pandas as pd
        from sklearn.preprocessing import StandardScaler

        self.scaler = StandardScaler()


        if self.set_type == 0:
            df_raw = pd.DataFrame(self.tr)
        elif self.set_type == 1:
            df_raw = pd.DataFrame(self.tr)
        else:
            df_raw = pd.DataFrame(self.ts)

        # Generate synthetic timestamps for columns (assume time flows along axis=1)
        sync = df_raw.copy().T
        base_time = pd.Timestamp('2025-01-01 00:00:00')

        timestamps = [base_time + pd.Timedelta(hours=i) for i in range(sync.shape[1])]
        sync.columns = timestamps  # assign synthetic timestamps to columns

        # Split sensor columns (rows) into train/val/test sets
        #sensor_indices = list(df_raw.index)
        #num_sensors = len(sensor_indices)
        #num_sensors_per_set = num_sensors // 3

        if self.set_type != 2:
            M, N = self.tr.shape
            split_idx = M // 2
            cont_tr = self.tr[:split_idx]  # first 1/3
            heal_tr = self.tr[split_idx:]  # remaining 2/3
            self.m_tr_cont = self.m_tr[:split_idx]  #
            self.m_tr_heal = self.m_tr[split_idx:]  #

            self.train_sensors = heal_tr
            self.val_sensors = cont_tr

            if self.set_type == 0:
                df_data = pd.DataFrame(self.train_sensors)
            elif self.set_type == 1:
                df_data = pd.DataFrame(self.val_sensors)
        else:
            M, N = self.ts.shape
            self.test_sensors = self.ts
            df_data = pd.DataFrame(self.test_sensors)

        # Apply standard scaling
        if self.scale:
            self.scaler.fit(df_data.values)
            data = self.scaler.transform(df_data.values)
        else:
            data = df_data.values

        # Create timestamp dataframe for time feature extraction
        df_stamp = pd.DataFrame({'date': sync.columns})
        df_stamp['date'] = pd.to_datetime(df_stamp['date'])

        # Encode time features
        if self.timeenc == 0:
            df_stamp['month'] = df_stamp['date'].apply(lambda row: row.month)
            df_stamp['day'] = df_stamp['date'].apply(lambda row: row.day)
            df_stamp['weekday'] = df_stamp['date'].apply(lambda row: row.weekday())
            df_stamp['hour'] = df_stamp['date'].apply(lambda row: row.hour)
            data_stamp = df_stamp.drop(['date'], axis=1).values
        elif self.timeenc == 1:
            data_stamp = time_features(pd.to_datetime(df_stamp['date'].values), freq=self.freq)
            data_stamp = data_stamp.transpose(1, 0)
        self.data_stamp = data_stamp

        # Store processed data
        self.data_x = data
        self.data_y = data

        if self.verbose:
            print(f"{self.data_x.shape = }, {self.data_y.shape = }, {self.data_stamp.shape = }")


    def __getitem__(self, index):

        s_begin = index
        s_end = s_begin + self.seq_len
        r_begin = s_end - self.label_len
        r_end = r_begin + self.label_len + self.pred_len

        seq_x = self.data_x[s_begin:s_end]
        seq_y = self.data_y[r_begin:r_end]
        seq_x_mark = self.data_stamp[s_begin:s_end]
        seq_y_mark = self.data_stamp[r_begin:r_end]

        # Get the appropriate mask
        if self.set_type == 0:  # train
            mask = self.m_tr_heal[s_begin:s_end]
        elif self.set_type == 1:  # val
            mask = self.m_tr_cont[s_begin:s_end]  # assuming same mask is used
        else:  # test
            mask = self.m_ts[s_begin:s_end]

        if self.verbose:
            print(f"Index {index} shapes: x={seq_x.shape}, y={seq_y.shape}, x_mark={seq_x_mark.shape}, y_mark={seq_y_mark.shape}, mask={mask.shape}, {self.seq_len = }, {self.patch_size = }, {self.pred_len = }, {self.label_len = }")

        return seq_x, seq_y, seq_x_mark, seq_y_mark, mask

    def __len__(self):
        return len(self.data_x) - self.seq_len - self.pred_len + 1

    def inverse_transform(self, data):
        return self.scaler.inverse_transform(data)