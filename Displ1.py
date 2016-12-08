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
    sact = Displ.fromfile("seact.log",N=88)
    wsu = Displ.fromfile("mowsu.log",N=88)
    dsu = Displ.fromfile("modsu.log",N=88)
    mact = Displ.fromfile("moact.log",N=88)
    plt.figure(fig)
    ax1 = plt.subplot(1,3,1)
    Displ.imshow(sact,axes=ax1,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    ax2 = plt.subplot(1,3,2)
    Displ.imshow(wsu,axes=ax2,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    ax3 = plt.subplot(1,3,3)
    Displ.imshow(mact,axes=ax3,r0 = r0,r1 = r1,u0 = u0,u1 = u1)
    plt.tight_layout()
    return sact,wsu,dsu,mact
