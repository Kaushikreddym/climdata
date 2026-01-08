# ===============================================================================================================
# SOURCE: https://github.com/caow13/BRITS
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://papers.nips.cc/paper_files/paper/2018/hash/734e6bfcd358e25ac1db0a4241b95651-Abstract.html
# ===============================================================================================================



import ujson as json
import numpy as np

import torch
from torch.utils.data import Dataset, DataLoader

class MySet(Dataset):
    def __init__(self, filename):
        super(MySet, self).__init__()
        self.content = open(filename).readlines()

        indices = np.arange(len(self.content))
        val_indices = np.random.choice(indices, len(self.content) // 5)

        self.val_indices = set(val_indices.tolist())

    def __len__(self):
        return len(self.content)

    def __getitem__(self, idx):
        rec = json.loads(self.content[idx])
        if idx in self.val_indices:
            rec['is_train'] = 0
        else:
            rec['is_train'] = 1
        return rec

def collate_fn(recs):
    forward = map(lambda x: x['forward'], recs)
    backward = map(lambda x: x['backward'], recs)

    def to_tensor_dict(recs):
        values = torch.FloatTensor(list(map(lambda r: r['values'], recs)))
        masks = torch.IntTensor(list(map(lambda r: r['masks'], recs)))
        deltas = torch.FloatTensor(list(map(lambda r: r['deltas'], recs)))
        forwards = torch.FloatTensor(list(map(lambda r: r['forwards'], recs)))
        evals = torch.FloatTensor(list(map(lambda r: r['evals'], recs)))
        eval_masks = torch.IntTensor(list(map(lambda r: r['eval_masks'], recs)))

        return {
            'values': values,
            'forwards': forwards,
            'masks': masks,
            'deltas': deltas,
            'evals': evals,
            'eval_masks': eval_masks
        }

    ret_dict = {
        'forward': to_tensor_dict(list(forward)),
        'backward': to_tensor_dict(list(backward))
    }

    ret_dict['labels'] = torch.FloatTensor(list(map(lambda x: x['label'], recs)))
    ret_dict['is_train'] = torch.FloatTensor(list(map(lambda x: x['is_train'], recs)))

    return ret_dict

def get_loader(filename, batch_size = 16, shuffle = False, num_workers=0):
    data_set = MySet(filename)

    data_iter = DataLoader(dataset = data_set, batch_size = batch_size, num_workers = num_workers, shuffle = shuffle, pin_memory = True, collate_fn = collate_fn)

    return data_iter
