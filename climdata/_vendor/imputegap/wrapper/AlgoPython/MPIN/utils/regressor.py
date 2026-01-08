# ===============================================================================================================
# SOURCE: https://github.com/XLI-2020/MPIN
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://www.vldb.org/pvldb/vol17/p345-li.pdf
# ===============================================================================================================


import torch
import torch.nn as nn
import numpy as np

def get_activation(activation):
    if activation == 'relu':
        return torch.nn.ReLU()
    elif activation == 'prelu':
        return torch.nn.PReLU()
    elif activation == 'tanh':
        return torch.nn.Tanh()
    elif activation == 'sigmoid':
        return torch.nn.Sigmoid()
    elif (activation is None) or (activation == 'none'):
        return torch.nn.Identity()
    else:
        raise NotImplementedError


class MLPNet(torch.nn.Module):
    def __init__(self,
         		input_dims, output_dim,
                 hidden_dim=32,
         		hidden_layer_sizes=(64,),
         		hidden_activation='relu',
         		output_activation=None,
                dropout=0.05):
        super(MLPNet, self).__init__()

        layers = nn.ModuleList()
        input_dim = np.sum(input_dims)
        # for layer_size in hidden_layer_sizes:
        # 	hidden_dim = layer_size
        # 	layer = nn.Sequential(
        # 				nn.Linear(input_dim, hidden_dim),
        # 				get_activation(hidden_activation),
        # 				nn.Dropout(dropout),
        # 				)
        # 	layers.append(layer)
        # 	input_dim = hidden_dim

        # layer = nn.Sequential(
        #                 nn.Linear(input_dim, hidden_dim),
        #     				get_activation(hidden_activation),
        #     				nn.Dropout(dropout),
        # 				nn.Linear(hidden_dim, output_dim),
        # 				get_activation(output_activation),
        # 				)

        layer = nn.Sequential(
            nn.Linear(input_dim, output_dim),
            get_activation(output_activation),
        )
       	layers.append(layer)
       	self.layers = layers

    def forward(self, inputs):
    	if torch.is_tensor(inputs):
    		inputs = [inputs]
    	input_var = torch.cat(inputs,-1)
    	for layer in self.layers:
    		input_var = layer(input_var)
    	return input_var




