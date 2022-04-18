import numpy as np
import matplotlib.pyplot as plt
import moosh as ms


def StructureBragg(periods):
    Eps = np.tile([2.25, 4.41], (1, periods))
    Eps = np.insert(Eps, 0, 1)
    Eps = np.append(Eps, 1)

    Mu = np.ones((2 * periods + 2))

    Height = np.tile([360 / (4 * 1.5), 360 / (4 * 2.1)], (1, periods))
    Height = np.insert(Height, 0, 1600)
    Height = np.append(Height, 100)

    return Eps, Mu, Height


struture = StructureBragg(10)

coef = ms.Reflection(structure=struture, polarisation=0)
coef.spectrum(0)
