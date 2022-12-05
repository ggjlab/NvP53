#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append("/public/home/guogjgroup/ggj/JiaqiLi/NvTK/")
print(sys.path)


# In[2]:


import h5py, os, argparse, logging, time

import numpy as np
import pandas as pd

import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

from NvTK import Trainer
from NvTK.Model.Publications import DeepSEA

from NvTK.Evaluator import calculate_roc, calculate_pr
from NvTK.Evaluator import show_auc_curve, show_pr_curve

from NvTK.Explainer import get_activate_W, meme_generate, save_activate_seqlets
from NvTK.Explainer import seq_logo, plot_seq_logo


# In[3]:


import seaborn as sns
import matplotlib.pyplot as plt


# In[4]:


os.makedirs("./Log", exist_ok=True)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=time.strftime('./Log/log_NvP53.%m%d.%H:%M:%S.txt'),
                    filemode='w')

# args
parser = argparse.ArgumentParser()
parser.add_argument("data_atac")
parser.add_argument("data_rna")
parser.add_argument("--gpu-device", dest="device_id", default="0")
args = parser.parse_args(['../dataset/Total-20220906/Dataset.pmat.pb.2kbp.20220921.h5', 
                          '../dataset/Total_RNA-20220906/Dataset.pb.10kbp.20220921.h5', 
                          '--gpu-device', '0'])
logging.info(args)


# In[5]:


## change device
os.environ["CUDA_VISIBLE_DEVICES"] = args.device_id
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# In[9]:

# unpack datasets
h5file = h5py.File(args.data_atac, 'r')

X_atac = h5file["X"][:].swapaxes(-1,1).astype(np.float32)
y_atac = h5file["y_orig"][:].astype(np.float32)

peaks = h5file["peak"][:].astype(str)
peaks_meta = h5file["peak_meta"][:].astype(str)

anno_atac = h5file["anno"][:].astype(str)
tasknames_atac = h5file["cluster_name"][:].astype(str)

h5file.close()

# In[10]:


t = 0.1
y_atac = np.where(y_atac > t, 1, 0).astype(np.float32)


# In[ ]:





# In[11]:


tasks_use = ['P53_1', 'P53_10', 'P53_11', 'P53_12', 'P53_13', 'P53_14',
       'P53_15', 'P53_16', 'P53_17', 'P53_18', 'P53_19', 'P53_2',
       'P53_20', 'P53_22', 'P53_24', 'P53_25', 'P53_3', 'P53_4', 'P53_5',
       'P53_6', 'P53_7', 'P53_8', 'P53_9', 'WT_1', 'WT_10', 'WT_11',
       'WT_12', 'WT_13', 'WT_15', 'WT_2', 'WT_20', 'WT_21', 'WT_22',
       'WT_23', 'WT_24', 'WT_25', 'WT_3', 'WT_4', 'WT_5', 'WT_6', 'WT_7',
       'WT_8', 'WT_9']

task_mask = [t in tasks_use for t in tasknames_atac]
tasknames_atac = tasknames_atac[task_mask]
y_atac = y_atac[:, task_mask]
logging.info(y_atac.shape)
logging.info(task_mask)


# In[12]:


# unpack anno
n_tasks_atac = y_atac.shape[-1]
n_length_atac = X_atac.shape[-1]

n_tasks_atac, tasknames_atac.shape, n_length_atac


# In[13]:


y_pos = np.sum(y_atac, axis=0)
y_pos_rate = y_pos / y_atac.shape[0]

plt.figure(figsize=(16,4))

plt.subplot(1, 2, 1)
sns.boxplot(y_pos_rate)
plt.xlabel("Positive Rate (y_pos / n_peaks)")
plt.ylabel("sample")

plt.subplot(1, 2, 2)
sns.barplot(x=tasknames_atac, y=y_pos_rate)
plt.xticks(rotation=90, size=8)
plt.xlabel("cluster")
plt.ylabel("Positive Rate (y_pos / n_peaks)")

plt.savefig("./boxplot.PosRate.y_used_atac.pdf")
plt.show()
plt.close()


# In[14]:


ko_mask_atac = ['P53' in n for n in tasknames_atac]
wt_mask_atac = ['WT' in n for n in tasknames_atac]
n_tasks_ko_atac=np.sum(ko_mask_atac); n_tasks_wt_atac=np.sum(wt_mask_atac)
n_tasks_ko_atac, n_tasks_wt_atac


# In[15]:


mask = peaks[:,-1]
train_idx, val_idx, test_idx = mask=='train', mask=='val', mask=='test'

x_train_atac, x_val_atac, x_test_atac = X_atac[train_idx], X_atac[val_idx], X_atac[test_idx]
y_train_atac, y_val_atac, y_test_atac = y_atac[train_idx], y_atac[val_idx], y_atac[test_idx]

np.sum(train_idx), np.sum(val_idx), np.sum(test_idx)


# In[ ]:





# In[16]:

# unpack datasets
h5file = h5py.File(args.data_rna, 'r')

X_rna = h5file["X"][:].swapaxes(-1,1).astype(np.float32)
y_rna = h5file["y"][:].astype(np.float32)

geneused = h5file["geneused"][:].astype(str)
geneused_meta = h5file["geneused_meta"][:].astype(str)

# anno = h5file["anno"][:].astype(str)
# anno_meta = h5file["anno_meta"][:].astype(str)
tasknames_rna = h5file["cluster_name"][:].astype(str)

h5file.close()

# In[17]:


# unpack anno
n_tasks_rna = y_rna.shape[-1]
n_length_rna = X_rna.shape[-1]

n_tasks_rna, tasknames_rna.shape, n_length_rna


# In[18]:


import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 4))
sns.heatmap(y_rna[10000:11000, :], cmap='Greys')
plt.show()
plt.close()


# In[19]:


y_rna = np.where(y_rna>1, 1., 0.)
y_rna = y_rna.astype(np.float32)
y_rna[:5,:10]


# In[20]:


y_pos = np.sum(y_rna, axis=0)
y_pos_rate = y_pos / y_rna.shape[0]

plt.figure(figsize=(16,4))

plt.subplot(1, 2, 1)
sns.boxplot(y_pos_rate)
plt.xlabel("Positive Rate (y_pos / n_genes)")
plt.ylabel("sample")

plt.subplot(1, 2, 2)
sns.barplot(x=tasknames_rna, y=y_pos_rate)
plt.xticks(rotation=90, size=8)
plt.xlabel("cluster")
plt.ylabel("Positive Rate (y_pos / n_genes)")

plt.savefig("./boxplot.PosRate.y_used_rna.pdf")
plt.show()
plt.close()


# In[21]:


ko_mask_rna = ['P53' in n for n in tasknames_rna]
wt_mask_rna = ['WT' in n for n in tasknames_rna]
n_tasks_ko_rna=np.sum(ko_mask_rna); n_tasks_wt_rna=np.sum(wt_mask_rna)
n_tasks_ko_rna, n_tasks_wt_rna


# In[22]:


mask = geneused[:,-1]
train_idx, val_idx, test_idx = mask=='train', mask=='val', mask=='test'

x_train_rna, x_val_rna, x_test_rna = X_rna[train_idx], X_rna[val_idx], X_rna[test_idx]
y_train_rna, y_val_rna, y_test_rna = y_rna[train_idx], y_rna[val_idx], y_rna[test_idx]

np.sum(train_idx), np.sum(val_idx), np.sum(test_idx)


# In[ ]:





# In[23]:


from torch.utils.data import Dataset, DataLoader

class NvDataset(Dataset):
    def __init__(self, X, xt, y_wt, y_ko):
        super().__init__()
        self.X = X
        self.xt = xt
        self.y_wt = y_wt
        self.y_ko = y_ko
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, index):
        x = self.X[index]
        xt = self.xt
        y_wt = self.y_wt[index]
        y_ko = self.y_ko[index]
        return x, xt, y_wt, y_ko


def collate_fn(data):
    X, Xt, Y_wt, Y_ko = [], [], [], []
    for x, xt, y_wt, y_ko in data:
        X.append(x)
        Xt.append(xt)
        Y_wt.append(y_wt)
        Y_ko.append(y_ko)
    X = torch.tensor(np.array(X))
    Xt = np.unique(Xt)
    Y_wt = torch.tensor(np.array(Y_wt))
    Y_ko = torch.tensor(np.array(Y_ko))
    return [X, Xt], [Y_wt, Y_ko]


# In[24]:


train_ds_atac = NvDataset(x_train_atac, "atac", y_train_atac[:, wt_mask_atac], y_train_atac[:, ko_mask_atac])
train_ds_rna = NvDataset(x_train_rna, "rna", y_train_rna[:, wt_mask_rna], y_train_rna[:, ko_mask_rna])
train_ds = train_ds_atac + train_ds_rna

val_ds_atac = NvDataset(x_val_atac, "atac", y_val_atac[:, wt_mask_atac], y_val_atac[:, ko_mask_atac])
val_ds_rna = NvDataset(x_val_rna, "rna", y_val_rna[:, wt_mask_rna], y_val_rna[:, ko_mask_rna])
val_ds = val_ds_atac + val_ds_rna

test_ds_atac = NvDataset(x_test_atac, "atac", y_test_atac[:, wt_mask_atac], y_test_atac[:, ko_mask_atac])
test_ds_rna = NvDataset(x_test_rna, "rna", y_test_rna[:, wt_mask_rna], y_test_rna[:, ko_mask_rna])
test_ds = test_ds_atac + test_ds_rna


# In[25]:


# define data loader
batch_size = 50

train_loader_rna = DataLoader(train_ds_rna, batch_size=batch_size, shuffle=True, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)
validate_loader_rna = DataLoader(val_ds_rna, batch_size=batch_size, shuffle=False, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)
test_loader_rna = DataLoader(test_ds_rna, batch_size=batch_size, shuffle=False, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)


# In[26]:


batch_size = 500

train_loader_atac = DataLoader(train_ds_atac, batch_size=batch_size, shuffle=True, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)
validate_loader_atac = DataLoader(val_ds_atac, batch_size=batch_size, shuffle=False, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)
test_loader_atac = DataLoader(test_ds_atac, batch_size=batch_size, shuffle=False, collate_fn=collate_fn, num_workers=0, drop_last=False, pin_memory=False)


# In[ ]:





# In[27]:


import torch
from torch import nn
from NvTK import BasicModule
from NvTK.Modules.Transformer import TransformerEncoder

# define modules
KoScaner = nn.Sequential(
            nn.Conv1d(4, 8, 19, 1, 9), 
            nn.ReLU(),
            nn.MaxPool1d(20),
            )

WtScaner = nn.Sequential(
            nn.Conv1d(4, 8, 19, 1, 9), 
            nn.ReLU(),
            nn.MaxPool1d(20),
            )

class Residual(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x, **kwargs):
        return self.fn(x, **kwargs) + x

SharedScaner = nn.Sequential(
            nn.Conv1d(4, 128, 11, 1, 5), nn.MaxPool1d(5),
            Residual(nn.Sequential(nn.Conv1d(128, 128, 9, 1, 4), nn.ReLU())),
            nn.Conv1d(128, 256, 7, 1, 3), nn.MaxPool1d(2),
            Residual(nn.Sequential(nn.Conv1d(256, 256, 7, 1, 3), nn.ReLU())),
            nn.Conv1d(256, 504, 5, 1, 2), nn.MaxPool1d(2)
            )


class SeqScaner(BasicModule):
    def __init__(self):
        super().__init__()
        self.KoScaner = KoScaner
        self.WtScaner = WtScaner
        self.SharedScaner = SharedScaner
        
        self.seq_pool_size = 20 # * 5
        self.SharedScaner_size = 504
        self.Scaner_size = 8
        self.seq_channel_size = 512 #1024

    def forward(self, x, t='Wt'):
        x0 = self.SharedScaner(x)
        if t == 'Ko':
            x1 = self.KoScaner(x)
        elif t == 'Wt':
            x1 = self.WtScaner(x)
        return [x0, x1]


class MultiPredictor(BasicModule):
    def __init__(self, n_size, n_tasks_ko=3, n_tasks_wt=2):
        super().__init__()
        self.linear_ko = nn.Linear(n_size, n_tasks_ko)
        self.linear_wt = nn.Linear(n_size, n_tasks_wt)
        
    def forward(self, x, t='Wt'):
        if t == 'Ko':
            x = self.linear_ko(x)
        elif t == 'Wt':
            x = self.linear_wt(x)
        return x


# define models
class NvATAC(BasicModule):
    def __init__(self, n_seqs=200, n_tasks_ko=3, n_tasks_wt=2):
        super().__init__()
        self.SeqScaner = SeqScaner()
        n_len = n_seqs // self.SeqScaner.seq_pool_size
        n_size = n_len * self.SeqScaner.seq_channel_size
        self.predictor = MultiPredictor(n_size, n_tasks_ko, n_tasks_wt)
        self.activate = nn.Sigmoid()

    def forward(self, x, t='Wt'):
        x = self.SeqScaner(x, t) # return a list
        x = torch.concat(x, dim=1)
        x = x.reshape(x.size(0), -1)
        x = self.predictor(x, t)
        x = self.activate(x)
        return x


class NvRNA(BasicModule):
    def __init__(self, n_seqs=2000, n_tasks_ko=3, n_tasks_wt=2):
        super().__init__()
        self.SeqScaner = SeqScaner()
        
        n_len = n_seqs // self.SeqScaner.seq_pool_size
        n_d = self.SeqScaner.SharedScaner_size
        self.rnn = nn.LSTM(input_size=n_d, hidden_size = n_d // 2 , 
                           num_layers=4, batch_first=True, bidirectional=True)
        
        n_size = n_len * (n_d + self.SeqScaner.Scaner_size)
        self.predictor = MultiPredictor(n_size, n_tasks_ko, n_tasks_wt)
        self.activate = nn.Sigmoid()

    def forward(self, x, t='Wt'):
        x0, x1 = self.SeqScaner(x, t) # bs, hidden, seq_len
        x0 = x0.swapaxes(1, -1) # bs, seq_len, embedding(batch_size,seq_length,input_size)
        x0, _ = self.rnn(x0) # bs, seq_len, hidden*2
        x0 = x0.swapaxes(1, -1) # bs, hidden*2, seq_len
        x = torch.concat([x0, x1], dim=1)
        x = x.reshape(x.size(0), -1)
        x = self.predictor(x, t)
        x = self.activate(x)
        return x


class NvP53(BasicModule):
    def __init__(self, n_seqs_rna=2000, n_tasks_ko_rna=3, n_tasks_wt_rna=2,
                    n_seqs_atac=200, n_tasks_ko_atac=3, n_tasks_wt_atac=2):
        super().__init__()
        self.NvRNA = NvRNA(n_seqs_rna, n_tasks_ko_rna, n_tasks_wt_rna)
        self.NvATAC = NvATAC(n_seqs_atac, n_tasks_ko_atac, n_tasks_wt_atac)
        # Sharing SeqScaner parameters
        self.NvRNA.SeqScaner = self.NvATAC.SeqScaner = SeqScaner()

    def forward(self, inp):
        x, xt = inp
        self.switch_xt(xt)
        if self.xt == 'rna':
            o1 = self.NvRNA(x, t='Wt')
            o2 = self.NvRNA(x, t='Ko')
        elif self.xt == 'atac':
            o1 = self.NvATAC(x, t='Wt')
            o2 = self.NvATAC(x, t='Ko')
        return [o1, o2, self.xt]

    def switch_xt(self, xt='atac'):
        self.xt = xt


# In[28]:


n_length_rna, n_tasks_ko_rna, n_tasks_wt_rna, n_length_atac, n_tasks_ko_atac, n_tasks_wt_atac


# In[29]:


model = NvP53(n_length_rna, n_tasks_ko_rna, n_tasks_wt_rna,
             n_length_atac, n_tasks_ko_atac, n_tasks_wt_atac)
model


# In[30]:


def criterion(outputs, targets):
    o1, o2, xt = outputs
    t1, t2 = targets
    loss = nn.BCELoss()
#     if xt == 'rna':
#         loss = nn.MSELoss()
#     else:
#         loss = nn.BCELoss()
    return loss(o1, t1) + loss(o2, t2)


# In[31]:


# x, y = next(iter(test_loader_rna))
# x = [i.to(device) if isinstance(i, torch.Tensor) else i for i in x]
# y = [i.to(device) if isinstance(i, torch.Tensor) else i for i in y]
# model.to(device)
# o = model(x)
# loss = criterion(o, y)
# loss.backward()


# In[32]:


optimizer = Adam(model.parameters(), lr=1e-6)
# criterion = nn.BCELoss().to(device)

trainer = Trainer(model, criterion, optimizer, device, patience=50,
                  use_tensorbord=False, evaluate_training=False)


# In[ ]:





# In[32]:


trainer.train_until_converge(train_loader_atac, validate_loader_atac, test_loader_atac, EPOCH=50)


# In[33]:


EPOCH = 500
for epoch in range(EPOCH):
    train_batch_loss_rna, train_loss_rna = trainer.train_per_epoch(train_loader_rna, epoch, verbose_step=20)
    val_batch_loss_rna, val_loss_rna, _, _ = trainer.predict(validate_loader_rna)
    test_batch_loss_rna, test_loss_rna, _, _ = trainer.predict(test_loader_rna)
    
    train_batch_loss_atac, train_loss_atac = trainer.train_per_epoch(train_loader_atac, epoch, verbose_step=20)
    val_batch_loss_atac, val_loss_atac, _, _ = trainer.predict(validate_loader_atac)
    test_batch_loss_atac, test_loss_atac, _, _ = trainer.predict(test_loader_atac)

    # logs
    _lr = trainer.optimizer.param_groups[0]['lr']

    train_metric = val_metric = test_metric = 0.5

    train_loss = train_loss_rna + train_loss_atac *100
    val_loss = val_loss_rna + val_loss_atac *100
    test_loss = test_loss_rna + test_loss_atac *100
    
    logging.info("Train\t Accuracy: ATAC Loss: %.4f\t\n" % (train_loss_rna))
    logging.info("Eval\t Accuracy: RNA Loss: %.4f\t\n" % (train_loss_atac))
    
    logging.info("Train\t Accuracy: %.4f\t Loss: %.4f\t\n" % (train_metric, train_loss))
    logging.info("Eval\t Accuracy: %.4f\t Loss: %.4f\t\n" % (val_metric, val_loss))
    
    trainer.logs['train_loss_list'].append(train_loss_rna)
    trainer.logs['val_loss_list'].append(val_loss_rna)
    trainer.logs['test_loss_list'].append(test_loss_rna)
    trainer.logs['train_metric_list'].append(train_loss_atac)
    trainer.logs['val_metric_list'].append(val_loss_atac)
    trainer.logs['test_metric_list'].append(test_loss_atac)
    trainer.logs['lrs'].append(_lr)
    trainer.show_trainer_log()

    # update best
    if val_loss < trainer.logs["best_val_loss"] or val_metric > trainer.logs["best_val_r"]:
        trainer.logs["best_val_loss"] = val_loss
        trainer.logs["best_val_r"] = val_metric
        trainer.logs["best_val_epoch"] = epoch

        logging.info("Eval\t Best Eval Accuracy: %.4f\t Loss: %.4f\t at Epoch: %d\t lr: %.8f\n" % (val_metric, val_loss, epoch, _lr))
        logging.info("Eval\t Test Accuracy: %.4f\t Loss: %.4f\n" % (test_metric, test_loss))
        
        trainer.save_checkpoint(best=True)

    # or checkpoint
    elif epoch % 20 == 1:
        trainer.save_checkpoint()

    # early stop
    if epoch >= trainer.logs["best_val_epoch"] + trainer.patience: #patience_epoch:
        break


# In[ ]:





# In[33]:


model = trainer.load_best_model("./Log/chekc_model.pth")


# In[34]:


# predict test-set
_, _, test_predictions_atac, test_targets_atac = trainer.predict(test_loader_atac)
test_predictions_atac = np.concatenate(test_predictions_atac, axis=-1)
test_targets_atac = np.concatenate(test_targets_atac, axis=-1)
test_targets_atac.shape, test_predictions_atac.shape


# In[35]:


test_predictions_atac.max()


# In[36]:


os.makedirs("./Test", exist_ok=True)

# metric test-set
fpr, tpr, roc_auc = calculate_roc(test_targets_atac, test_predictions_atac)
auroc = [roc_auc[k] for k in roc_auc.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

p, r, average_precision = calculate_pr(test_targets_atac, test_predictions_atac)
aupr = [average_precision[k] for k in average_precision.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

pd.DataFrame({"auroc":auroc, "aupr":aupr}, index=range(n_tasks_atac)).to_csv("Test/Metric.NvP53_atac.csv")


# In[37]:


np.mean(auroc), np.median(auroc)


# In[38]:


sel = range(500)#np.random.choice(range(test_targets_atac.shape[0]), 1000)

fig = plt.figure(figsize = (8, 8))

plt.subplot(2, 1, 1)
sns.heatmap(test_targets_atac[sel,], cmap='Greys')
plt.title("test_targets", fontsize=15)

plt.subplot(2, 1, 2)
sns.heatmap(test_predictions_atac[sel, ], cmap='Greys')
plt.title("test_predictions", fontsize=15)

fig.savefig("./Test/pt100.atac.png", format='png', dpi=300, bbox_inches='tight')
fig.show()


# In[39]:


show_auc_curve(fpr=fpr, tpr=tpr, roc_auc=roc_auc, save=True, fig_size=(5, 4))
show_auc_curve(fpr=fpr, tpr=tpr, roc_auc=roc_auc, save=False, fig_size=(5, 4))


# In[40]:


plt.figure(figsize=(6,4))
sns.boxplot(auroc)
plt.xlabel("auroc")
plt.savefig("./Test/boxplot.auroc_uq.pdf")
plt.show()
plt.close()


# In[41]:


plt.figure(figsize=(8,6))

plt.subplot(2, 1, 1)
sns.barplot(x=tasknames_atac, y=auroc)
plt.xticks(rotation=90, size=8)
plt.xlabel("cluster")
plt.ylabel("auroc")

plt.subplot(2, 1, 2)
sns.barplot(x=tasknames_atac, y=aupr)
plt.xticks(rotation=90, size=8)
plt.xlabel("cluster")
plt.ylabel("aupr")

# plt.savefig("./boxplot.PosRate.y.pb.0005.pdf")
plt.show()
plt.close()


# In[ ]:





# In[42]:


# predict test-set
_, _, test_predictions_rna, test_targets_rna = trainer.predict(test_loader_rna)
test_predictions_rna = np.concatenate(test_predictions_rna, axis=-1)
test_targets_rna = np.concatenate(test_targets_rna, axis=-1)
test_targets_rna.shape, test_predictions_rna.shape


# In[43]:


test_predictions_rna.max()


# In[44]:


genename = geneused[geneused[:,5]=="test", 4]
genename


# In[45]:


pd.DataFrame(test_targets_rna, index=genename, columns=tasknames_rna).to_csv("./Test/test_targets_rna.csv")
pd.DataFrame(test_predictions_rna, index=genename, columns=tasknames_rna).to_csv("./Test/test_predictions_rna.csv")


# In[46]:


os.makedirs("./Test", exist_ok=True)

# metric test-set
fpr, tpr, roc_auc = calculate_roc(test_targets_rna, test_predictions_rna)
auroc = [roc_auc[k] for k in roc_auc.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

p, r, average_precision = calculate_pr(test_targets_rna, test_predictions_rna)
aupr = [average_precision[k] for k in average_precision.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

pd.DataFrame({"auroc":auroc, "aupr":aupr}, index=range(n_tasks_rna)).to_csv("Test/Metric.NvP53_rna.csv")



# In[47]:


np.mean(auroc)


# In[48]:


sel = np.random.choice(range(test_targets_rna.shape[0]), 100)

fig = plt.figure(figsize = (8, 8))

plt.subplot(2, 1, 1)
sns.heatmap(test_targets_rna[sel,], cmap='Greys')
plt.title("test_targets", fontsize=15)

plt.subplot(2, 1, 2)
sns.heatmap(test_predictions_rna[sel, ], cmap='Greys')
plt.title("test_predictions", fontsize=15)

fig.savefig("./Test/pt100.rna.pdf", format='pdf', dpi=300, bbox_inches='tight')
fig.show()


# In[49]:


plt.figure(figsize=(6,4))
sns.boxplot(auroc)
plt.xlabel("pcc")
plt.savefig("./Test/boxplot.rna.pdf")
plt.show()
plt.close()


# In[50]:


def plot_tracks(tracks, interval, height=1.5, save=False):
    fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
    for ax, (title, y) in zip(axes, tracks.items()):
        ax.fill_between(np.linspace(interval.get("start"), interval.get("end"), num=len(y)), y)
        ax.set_title(title)
        sns.despine(top=True, right=True, bottom=True)
    ax.set_xlabel(str(interval))
    plt.tight_layout()
    plt.savefig(save) if save else None

# sns.heatmap(y_train_regulatory[:,:,0], cmap='Greys')


# In[51]:


best = pd.DataFrame({"pcc":auroc, "scc":aupr}, index=range(n_tasks_rna)).sort_values("pcc", ascending=False).head(5).index.values
best


# In[52]:


interval = dict(chrom="chr", start=0, end=12800000)

tracks = {}
for i in best:
    tracks['target_%d'%i] = test_targets_rna[:,i]
    tracks['pred_%d'%i] = test_predictions_rna[:,i]

plot_tracks(tracks, interval, height=1.5, save="Test/tracks_reg.pdf")


# In[ ]:





# ## Explain

# In[53]:


os.makedirs("./Explain", exist_ok=True)


# In[54]:


# hook
class ActivateFeaturesHook():
    def __init__(self, module):
        self.hook = module.register_forward_hook(self.hook_fn)
    def hook_fn(self, module, input, output):
        self.features = output.cpu().data.numpy()#.mean(-1)
    def get_features(self):
        return self.features
    def close(self):
        self.hook.remove()


def get_fmap(model, hook_module, data_loader, device=torch.device("cuda")):
    fmap, X = [], []
    model.eval()
    with torch.no_grad():
        activations = ActivateFeaturesHook(hook_module)
        for inp, _ in data_loader:
            x_tensor = [x.to(device) if isinstance(x, torch.Tensor) else x for x in inp]
            _ = model(x_tensor)
            X.append(x_tensor[0].cpu().numpy())
            fmap.append(activations.get_features())
        fmap = np.vstack(fmap)
        X = np.vstack(X)
        activations.close()
    return fmap, X


def get_activate_W_from_fmap(fmap, X, pool=1, threshold=0.99, motif_width=10, pad=0, axis=1):
    motif_nb = fmap.shape[1]
    X_dim, seq_len = X.shape[1], X.shape[-1]

    W=[]
    for filter_index in range(motif_nb):
        # find regions above threshold
        data_index, pos_index = np.where(fmap[:,filter_index,:] > np.max(fmap[:,filter_index,:], axis=axis, keepdims=True)*threshold)

        seq_align = []; count_matrix = []
        for i in range(len(pos_index)):
            # pad 1-nt
            start = pos_index[i] - 1 - pad
            end = start + motif_width + 2
            # handle boundary conditions
            if end > seq_len:
                end = seq_len
                start = end - motif_width - 2 
            if start < 0:
                start = 0 
                end = start + motif_width + 2 

            seq = X[data_index[i], :, start*pool:end*pool]
            seq_align.append(seq)
            count_matrix.append(np.sum(seq, axis=0, keepdims=True))

        seq_align = np.array(seq_align)
        count_matrix = np.array(count_matrix)

        # normalize counts
        seq_align = (np.sum(seq_align, axis=0)/np.sum(count_matrix, axis=0))*np.ones((X_dim, (motif_width+2)*pool))
        seq_align[np.isnan(seq_align)] = 0
        W.append(seq_align)

    W = np.array(W)
    return W


def get_activate_W(model, hook_module, data, pool=1, pad=0, threshold=0.99, motif_width=20, axis=1):
    fmap, X = get_fmap(model, hook_module, data)
    W = get_activate_W_from_fmap(fmap, X, pool, threshold, motif_width, pad=pad, axis=axis)
    return W


# In[55]:


# KoScaner
W = get_activate_W(model, model.NvATAC.SeqScaner.KoScaner[0], test_loader_atac, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.KoScaner.axis0.atac.npz", W=W)
meme_generate(W, output_file='./Explain/meme.KoScaner.axis0.atac.txt', prefix='Filter_')


# In[56]:


# WtScaner
W = get_activate_W(model, model.NvATAC.SeqScaner.WtScaner[0], test_loader_atac, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.WtScaner.atac.axis0.npz", W=W)
meme_generate(W, output_file='./Explain/meme.WtScaner.axis0.atac.txt', prefix='Filter_')


# In[57]:


# SharedScaner
W = get_activate_W(model, model.NvATAC.SeqScaner.SharedScaner[0], test_loader_atac, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.SharedScaner.axis0.atac.npz", W=W)
meme_generate(W, output_file='./Explain/meme.SharedScaner.axis0.atac.txt', prefix='Filter_')


# In[ ]:





# In[58]:


# KoScaner
W = get_activate_W(model, model.NvRNA.SeqScaner.KoScaner[0], test_loader_rna, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.KoScaner.axis0.rna.npz", W=W)
meme_generate(W, output_file='./Explain/meme.KoScaner.axis0.rna.txt', prefix='Filter_')


# In[59]:


# WtScaner
W = get_activate_W(model, model.NvRNA.SeqScaner.WtScaner[0], test_loader_rna, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.WtScaner.axis0.rna.npz", W=W)
meme_generate(W, output_file='./Explain/meme.WtScaner.axis0.rna.txt', prefix='Filter_')


# In[60]:


# SharedScaner
W = get_activate_W(model, model.NvRNA.SeqScaner.SharedScaner[0], test_loader_rna, threshold=0.99, motif_width=19, pad=9, axis=0)
np.savez_compressed("./Explain/W.SharedScaner.axis0.rna.npz", W=W)
meme_generate(W, output_file='./Explain/meme.SharedScaner.axis0.rna.txt', prefix='Filter_')


# In[ ]:





# ### calculate influence scores

# In[62]:


from NvTK.Explainer import deep_explain_layer_conductance
from NvTK.Explainer.Influence import channel_target_influence


# In[63]:


def foldchange(origin, modified):
    """caculate the fold change between modified and origin outputs."""
    return modified - origin


class ModifyOutputHook():
    def __init__(self, module):
        self.hook = module.register_forward_hook(self.hook_fn)
        self.channels = None
        self.channel = 0

    def hook_fn(self, module, input, output):
        for channel in self.channels:
            self.channel = channel
            if isinstance(module, torch.nn.modules.conv.Conv1d):
                output_channel = output[:,self.channel,:]
                output[:,self.channel,:] = torch.zeros_like(output_channel).to(output_channel.device)#output_channel.mean()
            elif isinstance(module, torch.nn.modules.linear.Linear):
                output_channel = output[:,self.channel]
                output[:,self.channel] = torch.zeros_like(output_channel).to(output_channel.device)#output_channel.mean()
            # logging.info(output_channel[:5].cpu().detach().numpy())
            # logging.info(output_channel.mean().cpu().detach().numpy())
        return output

    def step_channel(self, idx):
        if isinstance(idx, (list, tuple)):
            self.channels = idx
        elif isinstance(idx, int):
            self.channels = [idx]

    def get_current_channel(self):
        return self.channel

    def close(self):
        self.hook.remove()


def channel_target_influence(model, hook_module, data_loader, device=torch.device("cuda")):
    # criterion = torch.nn.BCELoss(reduction='none').to(device) # gene * cell
    target, pred_orig, loss_orig, pred_modified_foldchange = [], [], [], []

    # a normal feed-forward
    model.eval()
    with torch.no_grad():
        for inp, _ in data_loader:
            x_tensor = [x.to(device) if isinstance(x, torch.Tensor) else x for x in inp]
            output = model(x_tensor)[:2]
            # loss = criterion(output, t)

            # target.append(t.cpu().data.numpy())
            pred_orig.append(torch.concat(output, axis=-1).cpu().data.numpy())
            # loss_orig.append(loss.cpu().data.numpy())

        # target = np.vstack(target)
        pred_orig = np.vstack(pred_orig)
        # loss_orig = np.vstack(loss_orig)

        # feed-forward with ModifyOutputHook
        if isinstance(hook_module, torch.nn.modules.conv.Conv1d):
            out_channels = hook_module.out_channels # must hook on conv layer
        elif isinstance(hook_module, torch.nn.modules.linear.Linear):
            out_channels = hook_module.out_features # must hook on linear layer
        
        Modifier = ModifyOutputHook(hook_module)
        for idx in range(out_channels):
            logging.info("modifying channel_%d..." % idx)
            pred_modified, loss_modified = [], []
            Modifier.step_channel(idx)
            for inp, _ in data_loader:
                x_tensor = [x.to(device) if isinstance(x, torch.Tensor) else x for x in inp]
                output = model(x_tensor)[:2]
                
                pred_modified.append(torch.concat(output, axis=-1).cpu().data.numpy())
                # loss_modified.append(loss.cpu().data.numpy())
                
            pred_modified = np.vstack(pred_modified) 
            # loss_modified = np.vstack(loss_modified) 

            fc = foldchange(pred_orig, pred_modified).mean(0) # output_size
            # fc = foldchange(loss_orig, loss_modified).mean(0) # output_size
            pred_modified_foldchange.append(fc)

        Modifier.close()
    return np.vstack(pred_modified_foldchange)


# In[64]:

influe_conv = channel_target_influence(model, model.NvATAC.SeqScaner.SharedScaner[0], test_loader_atac)
influe_conv.shape

# In[65]:


np.savez_compressed("Explain/influe_diff.SharedScaner.atac.npz", data=influe_conv)


# In[ ]:





# In[66]:

influe_conv = channel_target_influence(model, model.NvRNA.SeqScaner.SharedScaner[0], test_loader_rna)
influe_conv.shape

# In[67]:


np.savez_compressed("Explain/influe_diff.SharedScaner.rna.npz", data=influe_conv)


# In[ ]:





# ## calculate layer conductance

# In[ ]:
imp = []
for tensor,_ in test_loader:
    tensor = tensor.to(device)
    conductance = deep_explain_layer_conductance(model, model.conv[0], tensor, n_tasks)
    imp.append(conductance)

# plot_label_neuron_importance(model, model.Embedding.conv, tensor, anno[:,0])

# In[ ]:


imp = np.hstack(imp)
imp.shape


# In[ ]:


np.savez_compressed("Explain/imp.npz", data=imp)


# In[ ]:


imp_mean = imp.mean(-1).mean(1)
imp_mean.shape


# In[ ]:


sns.clustermap(pd.DataFrame(imp_mean.T, columns=ann.pb_cluster_anno.values), # np.unique(anno_df.pb_cluster.values)), 
               cmap='vlag', standard_scale="row", row_cluster=True, figsize=(20, 5))


# In[ ]:


df = pd.DataFrame(imp.max(-1).mean(1)[:, :2].T, columns=ann.pb_cluster_anno.values)

plt.figure(figsize=(20, 3))
ax = sns.clustermap(df, cmap="Greys", figsize=(20, 5))
plt.savefig("./Explain/label_neuron_importance.ScanerP53.pdf")
plt.show()
plt.close()

