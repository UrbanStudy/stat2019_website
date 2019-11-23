''' Simple Classifier
Ethan Lew
10/10/19

A simple classifier designed to match the description of STAT 671 lecture 1.
'''
import numpy as np
from inspect import signature
from tqdm import tqdm

class LabeledData:
    ''' LabeledData
    LabeledData encapsulates the pairing between observations and {+/-1} labels
    '''
    def __init__(self):
        self._X = None
        self._Y = None

    @property
    def n(self):
        if self._X is None:
            return None
        else:
            return np.shape(self._X)[0]

    @property
    def d(self):
        if self._X is None:
            return None
        else:
            return np.shape(self._X)[1]

    def add_data(self, x, y):
        # Check that labels are valid
        unique = set(list(y))
        if unique != {-1, 1}:
            raise Exception("Y needs to contain just {-1, 1} binary labels!")

        # Add observations
        if self._X is None:
            self._X = x
            self._Y = y
        else:
            # Dimensionality Check
            if self.n != np.shape(x)[0]:
                raise Exception("New X data doesn't match dimensionality to the old data! (%s, %s)" %
                                (self.n, np.shape(x)[0]))
            self._X = np.vstack((self._X, x))
            self._Y = np.vstack((self._Y, y))

class PartitionData(LabeledData):
    '''PartitionedData
    PartitionedData stores a data frame and allows some of it to be used for training and validation
    '''
    def __init__(self):
        super(PartitionData, self).__init__()
        self._vX = None
        self._vY = None
        self._tX = None
        self._tY = None

    def partition(self, ratio):
        M = round(self.n * ratio)
        train = np.zeros((self.n), dtype=np.bool)
        train[0:M] = 1
        np.random.shuffle(train)
        self._tX = self._X[train]
        self._vX = self._X[~train]
        self._tY = self._Y[train]
        self._vY = self._Y[~train]

    @property
    def training(self):
        return self._tX, self._tY

    @property
    def validation(self):
        return self._vX, self._vY

class SimpleClassifier(LabeledData):
    ''' SimpleClassifier
    Encapsulates the classifier from class
    '''
    def __init__(self):
        super(SimpleClassifier, self).__init__()
        self._b = None
        self._xp = None
        self._xn = None
        self._mp = None
        self._mn = None
        self._k = np.dot

    @property
    def b(self):
        return self._b

    @property
    def k(self):
        return self._k

    @k.setter
    def k(self, func):
        if len(signature(func).parameters) != 2:
            raise Exception("kernel needs to have two arguments k(x, xp)!")
        else:
            self._k = func

    def train(self):
        self._xp = self._X[self._Y == 1, :]
        self._xn = self._X[self._Y == -1, :]
        self._mp = np.shape(self._xp)[0]
        self._mn = np.shape(self._xn)[0]
        b = 0
        total_calcs = np.shape(self._xp)[0]**2 + np.shape(self._xn)[0]**2
        with tqdm(total=total_calcs, desc="Training Classifier", bar_format="{l_bar}{bar} [ time left: {remaining} ]") as pbar:
            for xi in self._xp:
                for xj in self._xp:
                    b -= 1/(self._mp**2)*self._k(xi, xj)
                    pbar.update(1)
            for xi in self._xn:
                for xj in self._xn:
                    b += 1/(self._mn**2)*self._k(xi, xj)
                    pbar.update(1)
            b *= 0.5
            self._b = b

    def classify(self, x):
        g = self._b
        for xi in self._xp:
            g += 1/self._mp * self._k(x, xi)
        for xi in self._xn:
            g -= 1/self._mn * self._k(x, xi)
        return int(np.sign(g))


def risk(y, yp):
    return 1/(np.size(y))*np.sum(0.5*np.abs(y - yp))


