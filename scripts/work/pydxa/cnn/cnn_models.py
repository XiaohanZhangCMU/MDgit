import torch
import torch.nn as nn
from torch.autograd import Variable
import torch.nn.functional as F

class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        self.module_1 = nn.Sequential(
            nn.Conv2d(
                in_channels=3,
                out_channels=124,
                kernel_size=3,
                stride=1,
                padding=1,
            ),
            nn.BatchNorm2d(124),        
            nn.ReLU(),
#            nn.Conv2d(124,124,3,1,1),
#            nn.BatchNorm2d(124),        
#            nn.ReLU(),
            nn.Conv2d(124,124,3,1,1),
            nn.Dropout2d(p=0.1),
            nn.BatchNorm2d(124),        
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),                
        )
#        self.module_2 = nn.Sequential(
#                        nn.Conv2d(124, 248, 3, 1, 1),
#                        nn.BatchNorm2d(248),
#                        nn.ReLU(),
#                        nn.Conv2d(248,248, 3, 1, 1),
#                        nn.BatchNorm2d(248),
#                        nn.ReLU(),
#                        nn.Conv2d(248, 248, 3, 1, 1),
#                        nn.Dropout2d(p=0.1),
#                        nn.BatchNorm2d(248),
#                        nn.ReLU(),
#                        nn.MaxPool2d(2),
#                        )     


        self.out = nn.Linear(124, 2)

    def forward(self, x):
        x = self.module_1(x)
#        x = self.module_2(x)

        x = x.view(x.size(0), x.size(1), -1)
        x = x.mean(2)
        
        x = x.view(x.size(0),-1)
        output = self.out(x)
        output = F.log_softmax(output,dim=0)

        return output
    



