__author__ = 'Haohan Wang'

import numpy as np

def matchSNPs(X1, xname1, X2, xname2):
    ind1 = []
    ind2 = []

    xname1 = xname1.tolist()
    xname2 = xname2.tolist()

    newXname1 = []
    newXname2 = []

    for i in range(len(xname1)):
        if xname1[i] in xname2:
            ind1.append(i)
            i2 = xname2.index(xname1[i])
            ind2.append(i2)
            newXname1.append(xname1[i])
            newXname2.append(xname2[i2])

    Xnew1 = X1[:, np.array(ind1)]
    Xnew2 = X2[:, np.array(ind2)]

    return Xnew1, newXname1, Xnew2, newXname2