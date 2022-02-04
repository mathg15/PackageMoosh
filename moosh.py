import numpy as np
import matplotlib as plt


class Layer:

    def __init__(self, eps, mu, height):
        self.eps = eps
        self.mu = mu
        self.height = height

    def layer(self):
        return np.array([[self.eps, self.mu, self.height]])


class Structure:

    def __init__(self, *args):
        layer = args[0]
        for i in range(len(args) - 1):
            layer = np.append(layer, args[i + 1], axis=0)
        self.st = layer

    def structure(self):
        epsStruc = self.st[:, 0]
        muStruc = self.st[:, 1]
        heightStruc = self.st[:, 2]
        return epsStruc, muStruc, heightStruc

