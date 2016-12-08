import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time

import Displ

def displ1(r0 = 0,r1 = -1,u0 = 0,u1 = -1,fig = 1) :

    sact = None ; mwsu = None ; mdsu = None ; mact = None
    plt.figure(fig)
    if os.path.isfile("seact.log") :
        sact = Displ.fromfile("seact.log",N=88)
        ax1 = plt.subplot(1,3,1)
        Displ.imshow(sact,axes=ax1,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    if os.path.isfile("mowsu.log") :
        mwsu = Displ.fromfile("mowsu.log",N=88)
        ax2 = plt.subplot(1,3,2)
        Displ.imshow(mwsu,axes=ax2,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    if os.path.isfile("modsu.log") :
        mdsu = Displ.fromfile("modsu.log",N=88)
    if os.path.isfile("seact.log") :
        mact = Displ.fromfile("moact.log",N=88)
        ax3 = plt.subplot(1,3,3)
        Displ.imshow(mact,axes=ax3,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    plt.tight_layout()
    return sact,mwsu,mdsu,mact
