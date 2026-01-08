# ===============================================================================================================
# SOURCE: https://github.com/caow13/BRITS
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://papers.nips.cc/paper_files/paper/2018/hash/734e6bfcd358e25ac1db0a4241b95651-Abstract.html
# ===============================================================================================================

import json
import os

import numpy as np
import pandas as pd

def parse_delta(masks, dir_):
    if dir_ == "backward":
        masks = masks[::-1]

    deltas = []

    for t in range(len(masks)):
        if t == 0:
            deltas.append(np.ones(1))
        else:
            # if the value is not missing delta is 1 otherwise is 1 + prev delta
            deltas.append(np.ones(1) + (1 - masks[t]) * deltas[-1])

    return np.array(deltas)


def parse_rec(values, masks, evals, eval_masks, dir_):
    deltas = parse_delta(masks, dir_)

    forwards = pd.DataFrame(values).ffill().fillna(0.0).to_numpy()

    rec = {}

    rec["values"] = np.nan_to_num(values).tolist()
    rec["masks"] = masks.astype("int32").tolist()
    rec["evals"] = np.nan_to_num(evals).tolist()
    rec["eval_masks"] = eval_masks.astype("int32").tolist()
    rec["forwards"] = forwards.tolist()
    rec["deltas"] = deltas.tolist()

    return rec



def prepare_dat(input, output, m_train, m_test, m_val=None):

    os.makedirs(os.path.dirname(output), exist_ok=True)
    file = open(output, 'w')
    src = input

    for col in range(src.shape[0]):

        evals = []

        for v in src[col].tolist():
            evals.append([float(v)])

        evals = np.array(evals)
        shp = evals.shape

        evals = evals.reshape(-1)
        values = evals.copy()

        # Replace masks with external ones
        masks = m_train[col].reshape(-1)
        eval_masks = m_test[col].reshape(-1)

        masks = masks.astype(int)
        eval_masks = eval_masks.astype(int)

        evals = evals.reshape(shp)
        values = values.reshape(shp)
        masks = masks.reshape(shp)
        eval_masks = eval_masks.reshape(shp)

        rec = {"label": 0}

        rec["forward"] = parse_rec(values, masks, evals, eval_masks, dir_="forward")
        rec["backward"] = parse_rec(values[::-1], masks[::-1], evals[::-1], eval_masks[::-1], dir_="backward")

        rec = json.dumps(rec)

        file.write(rec + "\n")

    return rec
        
#end function
