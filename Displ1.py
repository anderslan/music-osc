import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time

import Displ

def mkmusicfile(inputasciifile,predata = None,row0 = 0,nrow = -1,col0 = 0,ncol = -1,nrepeat = 1,
                repeatspace = 0,lastspace = 0,fig = -1,outputasciifile = None) :
    NCOL = 88
    indata = np.loadtxt(inputasciifile)
    if nrow<0 : nrow = len(indata)

    predatalen = 0
    if predata!=None :
        predatalen = len(predata)

    totalrows = predatalen + nrepeat * nrow + (nrepeat-1) * repeatspace + lastspace
    utdata = np.zeros((totalrows,NCOL),np.float32)

    utdata[0:predatalen,:] = predata

    R = predatalen
    for rep in range(nrepeat) :
        for r in range(nrow) :
            utdata[R,col0:col0+ncol] = indata[r,col0:col0+ncol]
            R += 1
        if rep<nrepeat-1 :
            for r in range(repeatspace) : R += 1
        else :
            for r in range(lastspace) : R += 1

    if (fig>0) :
        plt.figure(fig)
        Displ.imshow(utdata,fig = fig)

    if outputasciifile!=None :
        np.savetxt(outputasciifile,utdata,fmt = "%1d")

    return utdata
