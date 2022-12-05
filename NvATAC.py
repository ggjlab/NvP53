import torch
from torch import nn
from NvTK import BasicModule
from NvTK.Modules import attend


class Residual(nn.Module):
    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def forward(self, x, **kwargs):
        return self.fn(x, **kwargs) + x

class NvAtac(BasicModule):
    def __init__(self, n_tasks):
        super().__init__()

        self.conv = nn.Sequential(
            nn.Conv1d(4, 128, 11, 1, 5), nn.LeakyReLU(0.2),
            Residual(nn.Sequential(nn.Conv1d(128, 128, 9, 1, 4), nn.LeakyReLU(0.2))),
            nn.MaxPool1d(4), nn.Dropout(0.2),

            nn.Conv1d(128, 256, 7, 1, 3), nn.LeakyReLU(0.2),
            Residual(nn.Sequential(nn.Conv1d(256, 256, 5, 1, 2), nn.LeakyReLU(0.2))),
            nn.MaxPool1d(4), nn.Dropout(0.2),

            nn.Conv1d(256, 512, 3, 1, 1), nn.LeakyReLU(0.2),
            Residual(nn.Sequential(nn.Conv1d(512, 512, 3, 1, 1), nn.LeakyReLU(0.2))),
            nn.MaxPool1d(4), nn.Dropout(0.2),

            nn.Conv1d(512, 1024, 3, 1, 1), nn.LeakyReLU(0.2),
            Residual(nn.Sequential(nn.Conv1d(1024, 1024, 3, 1, 1), nn.LeakyReLU(0.2))),
            nn.MaxPool1d(2), nn.Dropout(0.2),
            )
        
        self.linear = nn.Sequential(
            nn.Linear(1024, n_tasks+1),
            nn.LeakyReLU(0.2),
            nn.Linear(n_tasks+1, n_tasks),
            nn.Sigmoid()
            )

    def forward(self, x):
        x = self.conv(x) # bs, n_channels=1024, n_seq=4
        x = x.swapaxes(1, -1) # bs, n_seq=4, n_channels=1024
        x = self.linear(x)
        return x


class NvAtacPretrained(BasicModule):
    def __init__(self, n_tasks):
        super().__init__()
        self.NvAtac = NvAtac(1)
        self.linear = nn.Sequential(
            nn.Linear(1024, n_tasks+1),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(n_tasks+1, n_tasks),
            nn.Sigmoid()
            )

    def forward(self, x):
        x1 = self.conv(x)
        x2 = self.Scaner(x)
        x = torch.concat([x1, x2], dim=1)
        x = self.convs(x)
        x = x.reshape(x.size(0), -1)
        x = self.linear(x)
        return x

