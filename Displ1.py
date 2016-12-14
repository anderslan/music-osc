import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time
import datetime

import Displ

def displ1(r0 = 0,r1 = -1,u0 = 0,u1 = -1,fig = 1,mode = 'img') :

    slgi = None ; sact = None ; mwsu = None ; mdsu = None ; mact = None; mada = None;
    try:
        slgi_mtime = os.path.getmtime("selgi.log")
    except OSError:
        slgi_mtime = 0
    try:
        sact_mtime = os.path.getmtime("seact.log")
    except OSError:
        sact_mtime = 0
    try:
        mwsu_mtime = os.path.getmtime("mowsu.log")
    except OSError:
        mwsu_mtime = 0
    try:
        mdsu_mtime = os.path.getmtime("modsu.log")
    except OSError:
        mdsu_mtime = 0
    try:
        mact_mtime = os.path.getmtime("moact.log")
    except OSError:
        mact_mtime = 0
    try:
        mada_mtime = os.path.getmtime("moada.log")
    except OSError:
        mada_mtime = 0

    now = time.time()
    max_mtime = np.max([slgi_mtime,sact_mtime,mwsu_mtime,mdsu_mtime,mact_mtime,mada_mtime])

#    print max_mtime,slgi_mtime,sact_mtime,mwsu_mtime,mdsu_mtime,mact_mtime,mada_mtime
    if max_mtime <= now-60 : print "WARNING: All files are > 1 min old"

    nplt = 0
    if max_mtime - 1 <= slgi_mtime:
        slgi = Displ.fromfile("selgi.log",N=88)
        nplt += 1
    if max_mtime - 1 <= sact_mtime:
        sact = Displ.fromfile("seact.log",N=88)
        nplt += 1
    if max_mtime - 1 <= mwsu_mtime:
        mwsu = Displ.fromfile("mowsu.log",N=88)
        nplt += 1
    if max_mtime - 1 <= mdsu_mtime:
        mdsu = Displ.fromfile("modsu.log",N=88)
        nplt += 1
    if max_mtime - 1 <= mact_mtime:
        mact = Displ.fromfile("moact.log",N=88)
        nplt += 1
    if max_mtime - 1 <= mada_mtime:
        mada = Displ.fromfile("moada.log",N=88)
        nplt += 1

    plt.figure(fig)
    plt.clf()
    iplt = 0
    if slgi!=None:
        iplt += 1
        ax = plt.subplot(1,nplt,iplt)
        if mode=='img' :
            Displ.imshow(slgi,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(slgi,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
            ax.set_ylim(-0.1,1.1)
    if sact!=None:
        iplt += 1
        ax = plt.subplot(1,nplt,iplt)
        if mode=='img':
            Displ.imshow(sact,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(sact,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
            ax.set_ylim(-0.1,1.1)
    if mwsu!=None :
        iplt += 1
        ax = plt.subplot(1,3,2)
        if mode=='img' :
            Displ.imshow(mwsu,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(mwsu,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
    if mdsu!=None :
        iplt += 1
        ax = plt.subplot(1,nplt,iplt)
        if mode=='img' :
            Displ.imshow(mdsu,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(mdsu,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
    if mact!=None :
        iplt += 1
        ax = plt.subplot(1,nplt,iplt)
        if mode=='img' :
            Displ.imshow(mact,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(mact,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
            ax.set_ylim(-0.1,1.1)
    if mada!=None :
        iplt += 1
        ax = plt.subplot(1,nplt,iplt)
        if mode=='img' :
            Displ.imshow(mada,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
        elif mode=='plt' :
            Displ.plot(mada,axes=ax,r0 = r0,r1 = r1,u0 = u0,u1 = u1,clear=False)
    if (nplt>0) : plt.tight_layout()
    return slgi,sact,mwsu,mdsu,mact,mada
