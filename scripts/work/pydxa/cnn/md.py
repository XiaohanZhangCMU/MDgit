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
    train_list['cu']=['cu_0.npy']
    train_list['ni']=['ni_0.npy']
    train_list['si']=['si_0.npy']
    train_list['ge']=['ge_0.npy']

    test_list = {}
    test_list['cu']=['cu_0.npy']
    test_list['ni']=['ni_0.npy']
    test_list['si']=['si_0.npy']
    test_list['ge']=['ge_0.npy']

    labels = {}
    labels['cu']=1
    labels['ni']=1
    labels['si']=0
    labels['ge']=0

    def __init__(self, root, train=True):
        self.root = os.path.expanduser(root)
        self.train = train  # training set or test set
        self.total_train_data_points = 0
        self.total_test_data_points = 0
        # Let`s assume every data point has the same number of features, 3*base_features*base_features
        # If a data point has more than that, we throw an error
        # If a data point has less than that, we padd with zeros
        self.n_base_width = 120;
        self.n_base_height = 120;
        self.n_base_features = self.n_base_width*self.n_base_height*3
        # now load the picked numpy arrays
        if self.train:
            self.train_data = []
            self.train_labels = []
            for fentry in self.train_list:
                for file in fentry:
                    fo = open(os.path.join(root, file, 'rb')
                    entry = np.load(fo)
                    n_features = entry[2]
                    n_data_points = entry.size/(n_features+4) 
                    self.total_train_data_points += n_data_points
                    for i in range(n_data_points):
                        assert(n_features <= n_base_features), "Your data point {0} is too big! Truncate it!".format(i)
                        if (n_features < n_base_features):
                            zeros = [0]*(n_base_features-n_features)
                            self.train_data.append(list(entry[i*(n_features+4):(i+1)*(n_features+4)], zeros))
                        else:
                            self.train_data.append(entry[i*(n_features+4):(i+1)*(n_features+4)])
                        self.train_labels += self.labels[fentry.key]
                    fo.close()

            self.train_data = np.concatenate(self.train_data)
            self.train_data = self.train_data.reshape((self.total_data_points, 3, self.n_base_width, self.n_base_height))
        else:
            self.test_data = []
            self.test_labels = []
            for fentry in self.test_list:
                for file in fentry:
                    fo = open(os.path.join(root, file, 'rb')
                    entry = np.load(fo)
                    n_features = entry[2]
                    n_data_points = entry.size/(n_features+4) 
                    self.total_test_data_points += n_data_points
                    assert(n_features <= n_base_features), "One of your data point is too big! Truncate it!"
                    for i in range(n_data_points):
                        assert(n_features <= n_base_features), "Your data point {0} is too big! Truncate it!".format(i)
                        if (n_features < n_base_features):
                            zeros = [0]*(n_base_features-n_features)
                            self.test_data.append(list(entry[i*(n_features+4):(i+1)*(n_features+4)], zeros))
                        else:
                            self.test_data.append(entry[i*(n_features+4):(i+1)*(n_features+4)])
                        self.test_labels += self.labels[fentry.key]
                    fo.close()

            self.test_data = np.concatenate(self.test_data)
            self.test_data = self.test_data.reshape((self.total_test_data_points, 3, self.n_base_width, self.n_base_height))

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

