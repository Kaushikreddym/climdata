# ===============================================================================================================
# SOURCE: https://github.com/pbansal5/DeepMVI
#
# THIS CODE HAS BEEN MODIFIED TO ALIGN WITH THE REQUIREMENTS OF IMPUTEGAP (https://arxiv.org/abs/2503.15250),
#   WHILE STRIVING TO REMAIN AS FAITHFUL AS POSSIBLE TO THE ORIGINAL IMPLEMENTATION.
#
# FOR ADDITIONAL DETAILS, PLEASE REFER TO THE ORIGINAL PAPER:
# https://arxiv.org/abs/2103.01600
# ===============================================================================================================


from imputegap.wrapper.AlgoPython.DeepMVI.utils import *
from imputegap.wrapper.AlgoPython.DeepMVI.loader import *
from imputegap.wrapper.AlgoPython.DeepMVI.model import *

interval = 0

def train(model,train_loader,val_loader,device, max_epoch=1000, patience=2, lr=1e-3, seed=42, verbose=True):
    best_state_dict = model.state_dict()
    best_loss = float('inf')

    optim = torch.optim.Adam(model.parameters(),lr=lr)

    iteration = 0
    start_epoch = 0
    tolerance_epoch = 0
    train_error = 0

    for epoch in range(start_epoch,max_epoch):

        if verbose:
            print ("Starting Epoch : %d"%epoch)

        for inp_,mask,residuals,context_info in train_loader :
            inp_ = inp_.to(device).requires_grad_(True)
            loss = model(inp_,mask.to(device),residuals.to(device),context_info)
            optim.zero_grad()
            loss['mae'].backward()
            optim.step()
            iteration += 1
            train_error += float(loss['mae'].cpu())
            if (iteration % interval == 0):
                model.eval()
                loss_mre_num,count = 0,0
                with torch.no_grad():
                    for inp_,mask,residuals,context_info in val_loader :
                        loss = model.validate(inp_.to(device),mask.to(device), residuals.to(device),context_info)
                        loss_mre_num += (loss['loss_values']).sum()
                        count += len(loss['loss_values'])
                if (float(loss_mre_num)/count < 0.99*best_loss):
                    tolerance_epoch = 0
                    best_loss = float(loss_mre_num)/count
                elif (float(loss_mre_num)/count < best_loss):
                    best_state_dict = model.state_dict()
                    tolerance_epoch += 1
                else :
                    tolerance_epoch += 1
                if verbose:
                    print ('\tdone validation, Patience : ',tolerance_epoch)
                    print ('\tvalidation loss : ',float(loss_mre_num/count))
                    print ('\ttrain loss : ',float(train_error/interval))
                model.train()
                train_error = 0
                if (tolerance_epoch >= patience):
                    if verbose:
                        print ('\t\tEarly Stopping')
                    return best_state_dict
    return best_state_dict

def test(model,test_loader,val_feats,device):
    output_matrix = copy.deepcopy(val_feats)
    model.eval()
    with torch.no_grad():
        for inp_,mask,residuals,context_info in test_loader :
            loss = model.validate(inp_.to(device),mask.to(device),residuals.to(device),context_info,test=True)
            output_matrix[context_info[1][0]:context_info[1][0]+mask.shape[1],context_info[0][0,0]] = \
            np.where(mask.detach().cpu().numpy()[0],loss['values'].detach().cpu().numpy()[0],output_matrix[context_info[1][0]:context_info[1][0]+mask.shape[1],context_info[0][0,0]])
    model.train()
    return output_matrix


def transformer_recovery(input_feats, max_epoch=1000, patience=2, lr=1e-3, seed=42, verbose=True):
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.manual_seed(0)
    torch.cuda.manual_seed(0)
    torch.cuda.manual_seed_all(0)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    if verbose:
        print ('\nStart training DeepMVI\n')

    global interval

    mean = np.nanmean(input_feats,axis=0)
    std = np.nanstd(input_feats,axis=0)
    input_feats = (input_feats-mean)/std
    num_missing = 10*min(max(int(input_feats.shape[0]/100),1),500)

    train_feats, val_feats, val_points, test_points, block_size, kernel_size = make_validation(input_feats, num_missing=num_missing)
    
    if (block_size > 100):
        kernel_size = 20
    if block_size == 4:
        kernel_size = 4

    time_context = min(int(input_feats.shape[0]/2),30*kernel_size)

    use_embed= (not is_blackout(input_feats))
    use_context=(block_size <= kernel_size)
    use_local = (block_size < kernel_size)

    b_val = min(input_feats.shape[0], 16)
    batch_size = min(input_feats.shape[1] * int(input_feats.shape[0] / time_context), b_val)
    interval = 1000

    if verbose:
        print ('Block size is %d, kernel size is %d'%(block_size,kernel_size))
        print ('Use Kernel Regression : ',use_embed)
        print ('Use Context in Keys : ',use_context)
        print ('Use Local Attention : ',use_local, "\n")

        print(f"{train_feats.shape = }")
        print(f"{val_feats.shape = }")
        print(f"{val_points.shape = }")
        print(f"{test_points.shape = }")
        print(f"{batch_size = }\n")

    train_set = myDataset(train_feats,use_local,time_context = time_context)
    val_set = myValDataset(val_feats,val_points,False,use_local,time_context = time_context)
    test_set = myValDataset(val_feats,test_points,True,use_local,time_context = time_context)

    train_loader = torch.utils.data.DataLoader(train_set,batch_size = batch_size, drop_last = False,shuffle=False,collate_fn = my_collate)
    val_loader = torch.utils.data.DataLoader(val_set,batch_size = batch_size, drop_last = False,shuffle=False,collate_fn = my_collate)
    test_loader = torch.utils.data.DataLoader(test_set,batch_size = 1, drop_last = False,shuffle=False,collate_fn = my_collate)

    model = OurModel(sizes=[train_feats.shape[1]],kernel_size=kernel_size,block_size = block_size,nhead=2,time_len=train_feats.shape[0],use_embed=use_embed,use_context=use_context,use_local=use_local).to(device)
    model.std = torch.from_numpy(std).to(device)

    best_state_dict = train(model,train_loader,val_loader,device, max_epoch=max_epoch, patience=patience, lr=lr, seed=seed, verbose=verbose)
    model.load_state_dict(best_state_dict)

    matrix = test(model,test_loader,val_feats,device)
    matrix = (matrix*std)+mean
    return matrix
