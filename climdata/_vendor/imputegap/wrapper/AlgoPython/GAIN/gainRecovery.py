# coding=utf-8
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''Main function for UCI letter and spam datasets.
'''

# Necessary packages
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import warnings

import numpy as np

from imputegap.recovery.manager import TimeSeries
from imputegap.tools import utils
from imputegap.wrapper.AlgoPython.GAIN.gain import gain

warnings.simplefilter(action='ignore', category=FutureWarning)



def gainRecovery(miss_data_x, batch_size=-1, hint_rate=0.9, alpha=10, epoch=100, tr_ratio=0.9, verbose=True):
  '''Main function for UCI letter and spam datasets.
  
  Args:
    - data_name: letter or spam
    - miss_rate: probability of missing components
    - batch:size: batch size
    - hint_rate: hint rate
    - alpha: hyperparameter
    - iterations: iterations
    
  Returns:
    - imputed_data_x: imputed data
    - rmse: Root Mean Squared Error
  '''

  recov = np.copy(miss_data_x)
  m_mask = np.isnan(miss_data_x)
  input_data = np.copy(miss_data_x)

  if batch_size == -1:
      batch_size = utils.compute_batch_size(data=miss_data_x, min_size=4, max_size=32, verbose=verbose)

  if verbose:
      print(f"(IMPUTATION) GAIN\n\tMatrix: {input_data.shape[0]}, {input_data.shape[1]}\n\tbatch_size: {batch_size}\n\thint_rate: {hint_rate}\n\talpha: {alpha}\n\tepoch: {epoch}\n\ttr_ratio: {tr_ratio}\n")

  cont_data_matrix, mask_train, mask_test, mask_val, error = utils.dl_integration_transformation(input_data, tr_ratio=tr_ratio, inside_tr_cont_ratio=0.4, split_ts=1, split_val=0, nan_val=None, prevent_leak=False, offset=0.05, block_selection=False, seed=42, verbose=False)
  if error:
      return miss_data_x

  gain_parameters = {'batch_size': batch_size, 'hint_rate': hint_rate, 'alpha': alpha, 'iterations': epoch}

  max_trials = 10
  trial = 0
  tag = False

  while trial < max_trials:
      imputed_data_x = gain(cont_data_matrix, gain_parameters)

      if not np.all(np.isnan(imputed_data_x)):
          if tag:
              imputed_data_x = imputed_data_x.T
              break
          else:
              break

      cont_data_matrix = cont_data_matrix.T
      tag = ~tag
      trial += 1

  if verbose:
      print("All trials failed, returning last imputed result.")

  recov[m_mask] = imputed_data_x[m_mask]


  return recov  # Return last attempt even if NaN
