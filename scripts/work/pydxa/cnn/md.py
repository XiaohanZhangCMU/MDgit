from __future__ import print_function
import os
import os.path
import errno
import numpy as np
import sys
import torch
import torch.utils.data as data


class MD(data.Dataset):
    """ MD_ Dataset.
    Args:
        root (string): Root directory of dataset where all subdirectories of data files are stored
        train (bool, optional): If True, creates dataset from training set, otherwise
            creates from test set.
        transform (callable, optional): A function/transform that  takes in an atom-system
            and returns a transformed version. E.g, ``transforms.RandomCrop`` or padding
        target_transform (callable, optional): A function/transform that takes in the
            target and transforms it.
    """
    train_list = {}
    train_list['cu']=['data.npy']
    train_list['ni']=['data.npy']
    train_list['si']=['data.npy']
    train_list['ge']=['data.npy']

    test_list = {}
    test_list['cu']=['data.npy']
    test_list['ni']=['data.npy']
    test_list['si']=['data.npy']
    test_list['ge']=['data.npy']

    labels_1 = {}
    labels_1['cu']=1
    labels_1['ni']=1
    labels_1['si']=0
    labels_1['ge']=0

    def __init__(self, root, train=True):
        self.root = os.path.expanduser(root)
        self.train = train  # training set or test set
        self.total_train_data_points = 0
        self.total_test_data_points = 0
        # Let`s assume every data point has the same number of features, 3*base_features*base_features
        # If a data point has more than that, we throw an error
        # If a data point has less than that, we padd with zeros
        self.n_base_width = 128;
        self.n_base_height = 128;
        self.n_base_features = self.n_base_width*self.n_base_height*3

        train_start = 0
        train_end = 0.7

        test_start = 0.7
        test_end = 0.9
        tol = 1e-3; #zero tolerance for dislocation length

        # now load the picked numpy arrays
        if self.train:
            self.train_data = []
            self.train_labels = []
            for key, files in self.train_list.items():
                print(key)
                for npyfile in files:
                    print(npyfile)
                    fo = os.path.join(root, key, npyfile)
                    dps = np.load(fo)
                    n_features = int(dps[3])
                    print(n_features)
                    n_dps = dps.size/(n_features+4)
                    assert(len(dps.shape) == 1), "Shape of npydata has to be (N,)"
                    assert(n_features <= self.n_base_features), "Your data point {0} is too big! Truncate it!".format(n_features)
                    assert(n_dps-int(n_dps)==0), "Invalid npy file {0}".format(fo)
                    n_data_points = int(n_dps*(train_end-train_start))
                    print(n_data_points)
        
                    self.total_train_data_points += n_data_points
                    dps_validation = 0
                    for i in range(int(train_start*n_dps), int(train_end*n_dps)):
                        self.train_data.append(np.pad(dps[i*(n_features+4)+4:i*(n_features+4)+4+n_features], (0,self.n_base_features-n_features),'constant', constant_values=(0 )))  
#                        self.train_labels.append(self.labels_1[key])
                        labels_2 = 0
                        if np.abs(dps[i*(n_features+4)+2])>tol:
                            labels_2 = 1
                        self.train_labels.append(labels_2)
                        dps_validation += dps[i*(n_features+4)+3]
                    assert(dps_validation == n_data_points * n_features), "Number of features do not add up!"
            self.train_data = np.concatenate(self.train_data)
            self.train_data = self.train_data.reshape((self.total_train_data_points, 3, self.n_base_width, self.n_base_height))
            self.train_data = self.train_data.astype('float32')
#            self.train_data[:,:,:,:] = self.train_data
            print(self.train_data.shape)
            print(self.train_data.dtype)
        else:
            self.test_data = []
            self.test_labels = []
            for key, files in self.test_list.items():
                print(key)
                for npyfile in files:
                    print(npyfile)
                    fo = os.path.join(root, key, npyfile)
                    dps = np.load(fo)
                    n_features = int(dps[3])
                    print(n_features)
                    n_dps = dps.size/(n_features+4)
                    assert(len(dps.shape) == 1), "Shape of npydata has to be (N,)"
                    assert(n_features <= self.n_base_features), "Your data point {0} is too big! Truncate it!".format(n_features)
                    assert(n_dps-int(n_dps)==0), "Invalid npy file {0}".format(fo)
                    n_data_points = int(n_dps*(test_end-test_start))
                    print(n_data_points)
        
                    self.total_test_data_points += n_data_points
                    dps_validation = 0
                    for i in range(int(test_start*n_dps), int(test_end*n_dps)):
                        self.test_data.append(np.pad(dps[i*(n_features+4)+4:i*(n_features+4)+4+n_features], (0,self.n_base_features-n_features),'constant', constant_values=(0 )))  
#                        self.test_labels.append(self.labels_1[key])
                        labels_2 = 0
                        if np.abs(dps[i*(n_features+4)+2])>tol:
                            labels_2 = 1
                        self.test_labels.append(labels_2)
                        dps_validation += dps[i*(n_features+4)+3]
                    assert(dps_validation == n_data_points * n_features), "Number of features do not add up! {0} vs {1}".format(dps_validation, n_data_points * n_features)
            self.test_data = np.concatenate(self.test_data)
            self.test_data = self.test_data.reshape((self.total_test_data_points, 3,self.n_base_width, self.n_base_height))
            self.test_data = self.test_data.astype('float32')
            print(self.test_data.shape)
            print(self.test_data.dtype)
#            self.test_data[:,:,:,:] = self.test_data

    def __getitem__(self, index):
        """
        Args:
            index (int): Index
        Returns:
            tuple: (image, target) where target is index of the target class.
        """
        if self.train:
            img, target = self.train_data[index], self.train_labels[index]
        else:
            img, target = self.test_data[index], self.test_labels[index]

        img = torch.from_numpy(img) 
        return img, target

    def __len__(self):
        if self.train:
            return self.total_train_data_points
        else:
            return self.total_test_data_points

def normalize(x):
    mean = np.asarray([0.4914, 0.4822, 0.4465], np.float32)
    std = np.asarray([0.2023, 0.1994, 0.2010], np.float32)

    x = x.astype(np.float32)
    x = x / 255

    x -= mean.reshape((-1, 1, 1))
    x /= std.reshape((-1, 1, 1))

    return x

def hflip(x):
    if np.random.rand() > 0.5:
        x = x[:, :, ::-1]
    return x

def translate(x):
    new = np.zeros((3, 40, 40), np.float32)
    h = np.random.randint(9)
    w = np.random.randint(9)
    new[:, h:h + 32, w:w + 32] = x

    x = new[:, 4:36, 4:36]
    return x

