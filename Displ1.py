import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time

import Displ

def displ1(fig = 1) :
    lgi = Displ.fromfile("semolgi.log",N=88)
    wsu = Displ.fromfile("semowsu.log",N=88)
    act = Displ.fromfile("semoact.log",N=88)
    plt.figure(fig)
    ax1 = plt.subplot(1,3,1)
    Displ.imshow(lgi,axes=ax1)
    ax2 = plt.subplot(1,3,2)
    Displ.imshow(wsu,axes=ax2)
    ax3 = plt.subplot(1,3,3)
    Displ.imshow(act,axes=ax3)
    plt.tight_layout()
