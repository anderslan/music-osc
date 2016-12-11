import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time

import Displ

def mkdelay(dmin,dmax,n,tauzi) :
    if dmax<=dmin : 
        raise AssertionError("Mkdelay.mkdelay","Illegal dmax<=dmin");
    if n<=0: 
        raise AssertionError("Mkdelay.mkdelay","Illegal n<=0");
    ds = (dmax - dmin)/n
    tapsvec = []
    tauzivec = []
    for i in range(n) :
        tapsvec.append(dmin + (i+1)*ds)
        tauzivec.append(tauzi)
    return np.array(tapsvec),np.array(tauzivec)
