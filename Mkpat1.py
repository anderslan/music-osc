import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as distance
from numpy import ma
import select
import time

import Displ

def mkmusicdata(inputasciifile = None,predata = None,r0 = 0,r1 = -1,c0 = 0,c1 = -1,
                nrepeat = 1,repeatblank = 0,lastblank = 0,fig = -1,outasciifile = None) :
    NCOL = 88
    if inputasciifile!=None :
        indata = np.loadtxt(inputasciifile)
        if r1<0 : r1 = indata.shape[0]
        if c1<0 : c1 = indata.shape[1]

    predatalen = 0
    if predata!=None :
        predatalen = len(predata)

    totalrows = predatalen + nrepeat * (r1-r0) + (nrepeat-1) * repeatblank + lastblank
    utdata = np.zeros((totalrows,NCOL),np.float32)

    utdata[0:predatalen,:] = predata

    R = predatalen
    for rep in range(nrepeat) :
        for r in range(r0,r1) :
            utdata[R,c0:c1] = indata[r,c0:c1]
            R += 1
        if rep<nrepeat-1 :
            for r in range(repeatblank) : R += 1
        else :
            for r in range(lastblank) : R += 1

    if (fig>0) :
        plt.figure(fig)
        Displ.imshow(utdata,fig = fig)

    if outasciifile!=None :
        np.savetxt(outasciifile,utdata,fmt = "%1d")

    return utdata

def blankmusicdata(indata,r0 = 0,r1 = -1,u0 = 0,u1 = -1,outasciifile = None) :
    if r1<0 : r1 = len(indata)
    if u1<0 : u1 = indata.shape[1]
    
    data = np.array(indata)
    data[r0:r1,u0:u1] = 0.

    if outasciifile!=None :
        np.savetxt(outasciifile,data,fmt = "%1d")

    return data

def mergemusicdata(datalist,outasciifile = None) :
    alldata = datalist[0]
    for data in datalist[1:] :
        alldata = np.append(alldata,data,0)

    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdata1(outasciifile = None) :
    utdata1 = mkmusicdata("bwv772_100Hz.txt",r1 = 150,nrepeat = 5)
    utdata2 = blankmusicdata(utdata1,r0 = 400,r1 = 430,u1 = 44)
    alldata = mergemusicdata([utdata1,utdata2])
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdata2(outasciifile = None,u0 = -1) :
    alldata = mkmusicdata("bwv772_100Hz.txt",r1 = 100,nrepeat = 16)
    alldata = blankmusicdata(alldata,r0 = 400,r1=900)
    if u0<0 : u0 = alldata.shape[1]
    alldata = blankmusicdata(alldata,r0 = 900,u0=u0)
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdata3(outasciifile = None,u0 = -1,c1 = -1) :
    alldata1 = mkmusicdata("bwv772_100Hz.txt",r1 = 100,nrepeat = 14)
    alldata1 = blankmusicdata(alldata1,r0 = 400,r1=900)
    if u0<0 : u0 = alldata1.shape[1]
    alldata1 = blankmusicdata(alldata1,r0 = 900,u0=u0)
    alldata2 = mkmusicdata("bwv772_100Hz.txt",r0 = 100,r1=200,nrepeat = 12)
    alldata2 = blankmusicdata(alldata2,r1 = 300)
    alldata2 = blankmusicdata(alldata2,r0 = 700)
    if c1<0 : c1 = alldata1.shape[1]
    alldata3 = mkmusicdata("bwv772_100Hz.txt",r1 = 100,nrepeat = 5,c1=c1)
    alldata3 = blankmusicdata(alldata3,r0 = 300)
    alldata = mergemusicdata([alldata1,alldata2,alldata3])
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdata4(outasciifile = None) :
    alldata = mkmusicdata("bwv772_100Hz.txt",r1 = 800)
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkgliss(outasciifile = None) :
    data = np.zeros((2000,88),np.float32)
    for n in range(2000) :
        for i in range(11) :
            data[n,(n/10+8*i)%88] = 1
    if outasciifile!=None :
        np.savetxt(outasciifile,data,fmt = "%1d")
    return data
        
