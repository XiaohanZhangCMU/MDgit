'''
Factors that affect test accuracy:
1) Number of conv layer in each module.
2) Data augumentation. Implemented in cifar10.
3) Adaptive learning rate. 
4) Batch_size. 
5) Weigth_decay of SGD. 
6) Max pooling instead of fully connected layer.
7) Dropout layer.
8) Batch normalization.
9) Filter thickness.  
10) Adam vs SGD. 

To do: 
    write a dataloader for arbitrary data
    use cnn for dislocation detection problem
    problem w. utils.progressbar when submit jobs
'''

import sys
import numpy as np

import torch
import torchvision
import torch.nn as nn
from torch.autograd import Variable

import cnn_models
import cifar10
#from utils import progress_bar


N_EPOCH = 300
BATCH_SIZE = 128 #150
LR = 0.1 #Initial learning rate
DOWNLOAD_MNIST = False
DOWNLOAD_CIFAR10 = False
nnfile = 'cnn.pkl' 
nnparamfile = 'cnn.pkl.params'
use_cuda = torch.cuda.is_available()
def train_and_save( net, train_loader, test_loader, lr, N_EPOCH, nnfile, nnparamfile):
    loss_function = nn.CrossEntropyLoss()
    log_train = open('Log_Train_'+str(LR) +'_'+ str(BATCH_SIZE)+'.txt','a')
    log_valid = open('Log_Valid_'+str(LR) +'_'+ str(BATCH_SIZE)+'.txt','a')
    epoch_id = 0
   
    for epoch in range(N_EPOCH):
        train_loss = 0 
        total = 0
        correct = 0
        if (epoch_id < 1.0/3 *N_EPOCH):
            lr = LR
        elif (epoch_id < 2.0/3*N_EPOCH):
            lr = 0.01
        else:
            lr = 0.0028

        # train the current epoch
        # optimizer = torch.optim.Adam(net.parameters(), lr=lr)
        optimizer = torch.optim.SGD(net.parameters(), lr=lr, momentum=0.9, weight_decay=1e-4)
        for batch_idx, (x,y) in enumerate(train_loader):
            if (use_cuda): 
                x, y = x.cuda(), y.cuda()

            b_x, b_y = Variable(x), Variable(y)
            prediction = net(b_x)
            loss = loss_function(prediction, b_y)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            train_loss += loss.data[0]
            _, predicted = torch.max(prediction.data,1)
            total += b_y.size(0)
            correct += predicted.eq(b_y.data).cpu().sum()

        buff = 'epoch =' + str(epoch)+ ': train_loss: ' + str(train_loss/(batch_idx+1)) + ': train accuracy: ' + str(100.*correct/total)+'\n'
        log_train.write(buff)
        print(buff) 

        # test the current epoch
        net.eval() # switch net to 'test' mode
        test_loss = 0
        correct = 0
        total = 0
        loss_function = nn.CrossEntropyLoss()
        for batch_idx, (inputs, targets) in enumerate(test_loader):
            if use_cuda:
                inputs, targets = inputs.cuda(), targets.cuda()
            inputs, targets = Variable(inputs, volatile=True), Variable(targets)
            outputs = net(inputs)
            loss = loss_function(outputs, targets)
            test_loss += loss.data[0]
            _, predicted = torch.max(outputs.data, 1)
            total += targets.size(0)
            correct += predicted.eq(targets.data).cpu().sum()
            
        buff = 'epoch =' + str(epoch)+ ': test accuracy: ' + str(100.*correct/total)+'\n'
        print(buff)
       
        net.train() # switch net to 'train' mode
        epoch_id += 1

#   end of training
    torch.save(net, nnfile+'.'+str(LR)+'.'+str(BATCH_SIZE))
    torch.save(net.state_dict(), nnparamfile+'.'+str(LR)+'.'+str(BATCH_SIZE))
    log_train.close()
    log_valid.close()
    return net

# train to recognize MNIST data
def main_1(argv):
    train_data = torchvision.datasets.MNIST(
        root = '../mnist',
        train = True,
        transform = torchvision.transforms.ToTensor(),
        download=DOWNLOAD_MNIST
    )
    train_loader = torch.utils.data.DataLoader(dataset = train_data, batch_size=BATCH_SIZE,
                                   shuffle = True, num_workers = 2)
    test_data = torchvision.datasets.MNIST(root='../mnist/', train=False)
    validation_x = Variable(torch.unsqueeze(test_data.test_data, dim=1),volatile=True).type(torch.FloatTensor)[:2000]/255.
    validation_y = test_data.test_labels[:2000]

    net = cnn_models.CNN( )
    net = train_and_save( net, train_loader, validation_x, validation_y, LR, N_EPOCH, nnfile, nnparamfile)

# train to reconginize CIFAR10 data
def main(argv):
    train_data = cifar10.CIFAR10(
        root = '../cifar10',
        train = True,
        download=DOWNLOAD_CIFAR10
    )
    train_loader = torch.utils.data.DataLoader(dataset = train_data, batch_size=BATCH_SIZE,
                                               shuffle = True, num_workers = 2)
    test_data = cifar10.CIFAR10(root='../cifar10/', train=False, download=False)
    test_loader = torch.utils.data.DataLoader(dataset = test_data, batch_size = 100, shuffle=True, num_workers=2) 

    net = cnn_models.CNN( )
    if use_cuda:
        net.cuda()
        net = torch.nn.DataParallel(net, device_ids=range(torch.cuda.device_count())) 
        torch.backends.cudnn.enabled=True
        
    net = train_and_save( net, train_loader, test_loader, LR, N_EPOCH, nnfile, nnparamfile)

if __name__ == "__main__":
    main(sys.argv[1:])
