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
                out_channels=128,
                kernel_size=3,
                stride=1,
                padding=1,
            ),
            nn.BatchNorm2d(128),        
            nn.ReLU(),
            nn.Conv2d(128,128,3,1,1),
            nn.BatchNorm2d(128),        
            nn.ReLU(),
            nn.Conv2d(128,128,3,1,1),
            nn.Dropout2d(p=0.4),
            nn.BatchNorm2d(128),        
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),                
        )
        self.module_2 = nn.Sequential(
                        nn.Conv2d(128, 256, 3, 1, 1),
                        nn.BatchNorm2d(256),
                        nn.ReLU(),
                        nn.Conv2d(256,256, 3, 1, 1),
                        nn.BatchNorm2d(256),
                        nn.ReLU(),
                        nn.Conv2d(256, 256, 3, 1, 1),
                        nn.Dropout2d(p=0.4),
                        nn.BatchNorm2d(256),
                        nn.ReLU(),
                        nn.MaxPool2d(2),
                        )     

#        self.module_3 = nn.Sequential(
#            nn.Conv2d(256, 512, 3, 1, 1),
#            nn.BatchNorm2d(512),        
#            nn.ReLU(),
#            nn.Conv2d(512, 512, 3, 1, 1),
#            nn.BatchNorm2d(512),        
#            nn.ReLU(),
#            nn.Conv2d(512, 512, 3, 1, 1),
#            nn.Dropout2d(p = 0.4),
#            nn.BatchNorm2d(512),        
#            nn.ReLU(),
#            nn.MaxPool2d(2),                                
#        )
#
#        self.module_4 = nn.Sequential(
#            nn.Conv2d(512, 1024, 3, 1, 1),
#            nn.BatchNorm2d(1024),        
#            nn.ReLU(),
#            nn.Conv2d(1024, 1024, 3, 1, 1),
#            nn.BatchNorm2d(1024),        
#            nn.ReLU(),
#            nn.Conv2d(1024, 1024, 3, 1, 1),
#            nn.Dropout2d(p = 0.4),
#            nn.BatchNorm2d(1024),        
#            nn.ReLU(),
#            nn.MaxPool2d(2),                                
#        )

        #self.module_5 = nn.Sequential(
        #    nn.Conv2d(1024, 2048, 3, 1, 1),
        #    nn.BatchNorm2d(2048),        
        #    nn.ReLU(),
        #    nn.Conv2d(2048, 2048, 3, 1, 1),
        #    nn.BatchNorm2d(2048),        
        #    nn.ReLU(),
        #    nn.Conv2d(2048, 2048, 3, 1, 1),
#       #     nn.Dropout2d(p = 0.1),
        #    nn.BatchNorm2d(2048),        
        #    nn.ReLU(),
        #    nn.MaxPool2d(2),                                
        #)
        self.out = nn.Linear(256, 2)

    def forward(self, x):
        x = self.module_1(x)
        x = self.module_2(x)
#        x = self.module_3(x)
#        x = self.module_4(x)
#       x = self.module_5(x)

        x = x.view(x.size(0), x.size(1), -1)
        x = x.mean(2)
        
        x = x.view(x.size(0),-1)
        output = self.out(x)
        output = F.log_softmax(output,dim=0)

        return output
    
    torch.Size([2, 512])


