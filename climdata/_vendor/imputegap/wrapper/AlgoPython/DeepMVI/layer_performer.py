# ===============================================================================================================
# SOURCE: https://github.com/pbansal5/DeepMVI
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/abs/2103.01600
# ===============================================================================================================


from typing import Optional
import torch
from torch import Tensor


class FastCausalMultiheadAttention(torch.nn.Module):
    def __init__(self,d_query:int,d_key:int,d_value:int, nhead:int, backward=False):
        super(FastCausalMultiheadAttention, self).__init__()
        self.head_dim = d_value
        self.nhead = nhead
        self.scaling = float(d_query) ** -0.5
        self.d_query = d_query
        self.d_key = d_key
        self.d_value = d_value
        self.query_linear = torch.nn.Linear(d_query,nhead*d_query)
        self.key_linear = torch.nn.Linear(d_key,nhead*d_key)
        self.value_linear = torch.nn.Linear(d_value,nhead*d_value)
        self.out_linear = torch.nn.Linear(nhead*d_value,d_value)

        nb_features = int(d_query * math.log(d_query))
        self.projection_matrix = gaussian_orthogonal_random_matrix(nb_rows = nb_features, nb_columns = d_query, scaling = 0)
        self.kernel_fn = torch.nn.ReLU()
        self.causal_linear_fn = causal_linear_attention
        self.mask_size = 10
        self.backward = backward
        
    def forward(self,query: Tensor,key: Tensor, value: Tensor,attn_mask: Optional[Tensor] = None) -> Tensor:
        tgt_len, bsz, _ = query.size()
        query = self.query_linear(query.transpose(0,1)).view(bsz,tgt_len,self.nhead,self.d_query).transpose(1,2)
        key = self.key_linear(key.transpose(0,1)).view(bsz,tgt_len,self.nhead,self.d_key).transpose(1,2)
        value = self.value_linear(value.transpose(0,1)).view(bsz,tgt_len,self.nhead,self.d_value).transpose(1,2)
        
        query *= self.scaling
        
        device = query.device
        
        query = softmax_kernel(query, is_query = True, projection_matrix = self.projection_matrix, device = device)
        key = softmax_kernel(key, is_query = False, projection_matrix = self.projection_matrix, device = device)

        if (self.training):
            mask_len = torch.randint(low=1,high=self.mask_size+1,size=(1,))[0]
        else :
            mask_len = 1
        dummy = torch.zeros((bsz,self.nhead,mask_len,self.d_value)).to(query.device)
        if (not self.backward):
            temp = self.causal_linear_fn(query[:,:,mask_len:,:], key[:,:,:-mask_len,:], value[:,:,:-mask_len,:])
            temp = torch.cat([dummy,temp],dim=2)
        else :
            temp = self.causal_linear_fn(query[:,:,:-mask_len,:], key[:,:,mask_len:,:], value[:,:,mask_len:,:])
            temp = torch.cat([temp,dummy],dim=2)
        # if (self.training):
        #     mask_len = torch.randint(low=1,high=self.mask_size+1,size=(1,))[0]
        #     dummy = torch.zeros((bsz,self.nhead,mask_len,self.d_value)).to(query.device)
        #     if (not self.backward):
        #         temp = self.causal_linear_fn(query[:,:,mask_len:,:], key[:,:,:-mask_len,:], value[:,:,:-mask_len,:])
        #         temp = torch.cat([dummy,temp],dim=2)
        #     else :
        #         temp = self.causal_linear_fn(query[:,:,:-mask_len,:], key[:,:,mask_len:,:], value[:,:,mask_len:,:])
        #         temp = torch.cat([temp,dummy],dim=2)
        # else :
        #     temp = self.causal_linear_fn(query, key, value)
            
        value = self.out_linear(temp.transpose(1,2).contiguous().view(bsz,tgt_len,-1)).transpose(0,1)

        return value



class PerformerEncoderLayer(torch.nn.Module):
    def __init__(self,d_query,d_key,d_value,nhead, dim_feedforward=2048, dropout=0.1, backward=False,activation="relu"):
        super(PerformerEncoderLayer, self).__init__()
        self.backward = backward
        self.self_attn1 = FastCausalMultiheadAttention(d_query,d_key,d_value, nhead,backward=self.backward)
        # self.self_attn2 = FastCausalMultiheadAttention(d_query,d_key,d_value, nhead,backward=self.backward)
        self.self_attn2 = FastCausalMultiheadAttention(d_value,d_value,d_value, nhead,backward=self.backward)
        self.linear1 = torch.nn.Linear(d_value, dim_feedforward)
        self.dropout = torch.nn.Dropout(dropout)
        self.linear2 = torch.nn.Linear(dim_feedforward, d_value)
        self.norm = torch.nn.LayerNorm(d_value)
        self.activation = torch.nn.ReLU()
        
    def forward(self, query: Tensor, key: Tensor, value: Tensor,src_mask: Optional[Tensor] = None) -> Tensor:
        src = self.self_attn1(query, key, value, src_mask)
        src = self.linear2(self.dropout(self.activation(self.linear1(self.activation(src)))))
        src[torch.isnan(src)] = 0

        # src = self.self_attn2(query, key, src, src_mask)
        # src = self.linear2(self.dropout(self.activation(self.linear1(self.activation(src)))))
        # src[torch.isnan(src)] = 0

        src = self.self_attn2(src, src, src, src_mask)
        src = self.linear2(self.dropout(self.activation(self.linear1(self.activation(src)))))
        src[torch.isnan(src)] = 0

        return src

