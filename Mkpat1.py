import sys, os
import math
import numpy as np
import matplotlib.pyplot as plt
from random import randint
from numpy.random import poisson
import scipy.spatial.distance as distance
import select
import time

import Displ

NCOL = 88

def mkrandtimedpat_(H,U,udur,adur) :
    if udur<=0 :
        raise AssertionError("Timedpat.mkrandtimedpat","Illegal udur<=0")
    if H>1 :
        timedpats = mkrandtimedpat_(1,U,udur,adur)
        for h in range(1,H) :
            timedpat = mkrandtimedpat_(1,U,udur,adur)
            timedpats = np.concatenate([timedpats,timedpat],1)
        return timedpats
    itime = 0
    timedpat = np.zeros((adur,U),np.float32)
    while itime<adur :
        u = randint(0,U-1)
        dur = poisson(udur)
        for i in range(0,dur) :
            if adur<=itime : break
            timedpat[itime,u] = 1
            itime += 1
    return timedpat

def mkrandtimedpat(H,U,udur,adur,outasciifile = None) :
    pp = mkrandtimedpat_(H,U,udur,adur)
    if outasciifile!=None :
        np.savetxt(outasciifile,pp,fmt = "%1d")
    return pp

def mkmusicdata(indata,predata = None,r0 = 0,r1 = -1,c0 = 0,c1 = -1,
                nrepeat = 1,repeatblank = 0,lastblank = 0,fig = -1,outasciifile = None) :
    if type(indata)==str :
        indata = np.loadtxt(indata)
        if r1<0 : r1 = indata.shape[0]
        if c1<0 : c1 = indata.shape[1]
    else :
        if r1<0 : r1 = indata.shape[0]

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

def blankdata(nrow) :
    return np.zeros((nrow,NCOL))

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

def mkdata5(outasciifile = None) :
    indata = mkrandtimedpat(2,44,3,100)
    alldata = mkmusicdata(indata,r1 = 100,nrepeat = 14)
    alldata = blankmusicdata(alldata,r0 = 400)
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdata6(outasciifile = None,n = 50,dur = 2,partdur = 100,u0 = 50,nrepeat = 25,
            blankr0 = 400) :
    alldata = np.zeros((n*dur,NCOL))
    for p in range(n*dur) :
        alldata[p,22 + p/dur] = 1
        alldata[p,33 + p/dur] = 1
    alldata = mkmusicdata(alldata,nrepeat = nrepeat)
    alldata = blankmusicdata(alldata,r0 = blankr0,r1 = 250)
    alldata = blankmusicdata(alldata,r0 = 1000,r1 = 1500)
    alldata = blankmusicdata(alldata,r0 = 1100,r1 = 1900 ,u0 = u0)
    alldata = blankmusicdata(alldata,r0 = 1900)
    if outasciifile!=None :
        np.savetxt(outasciifile,alldata,fmt = "%1d")
    return alldata

def mkdashedstripes(data,r0 = 0, r1 = -1,dur = 1,gap = 0,uu = []) :
    if r1<0 : r1 = len(data)
    if len(data)<=r0 : raise AssertionError("Timedpat.mkdata7","Illegal len(data)<=r0")
    if len(data)<r1 : raise AssertionError("Timedpat.mkdata7","Illegal len(data)<r1")
    if r0>r1 : raise AssertionError("Timedpat.mkdata7","Illegal r0>r1)")
    for r in range(r0,r1) :
        for u in uu :
            if u<0 or u>len(data)-1 : raise AssertionError("Timedpat.mkdata7","Illegal u>len(data)")
            rn = float(r - r0)/(dur + gap - 1)
            data[r,u] = int(rn)%2
#            print r,rn,data[r,u]

def mkdata00(nrow,nrepeat = 5,repeatblank= 20,lastblank = 500,savetofile = False) :
    alldata = mkmusicdata("bwv772_100Hz.txt",r1 = nrow,nrepeat = nrepeat,repeatblank = repeatblank,
                          lastblank = lastblank)
    if savetofile :
        np.savetxt("data00.txt",alldata,fmt = "%1d")
    return alldata

def mkdata01(nrow = 200,nrepeat = 5,repeatblank= 20,lastblank = 400,c0 = 0,c1 = 40,savetofile = False) :
    data1 = mkmusicdata("bwv772_100Hz.txt",r1 = nrow,nrepeat = nrepeat,repeatblank = repeatblank,
                           lastblank = 500)
    data2 = mkmusicdata("bwv772_100Hz.txt",r1 = nrow,nrepeat = 2,lastblank = lastblank,c0 = c0,c1 = c1)
    alldata = mergemusicdata([data1,data2])
    if savetofile :
        np.savetxt("data01.txt",alldata,fmt = "%1d")
    return alldata

def mkdata02(savetofile = False) :
    alldata = mkdata6(n = 15,dur = 7,u0 = 35)
    if savetofile :
        np.savetxt("data02.txt",alldata,fmt = "%1d")
    return alldata

def mkdata03(nrow = 1000,savetofile = False) :
    alldata = np.zeros((nrow,NCOL))
    mkdashedstripes(alldata,r0 = 30,r1 = 210,uu = [0],dur = 25,gap = 5)
    mkdashedstripes(alldata,r0 = 500,r1 = 650,uu = [0],dur = 25,gap = 5)
#    mkdashedstripes(alldata,r0 = 30,r1 = nrow-30,uu = [20,22,24,40,42,44],dur = 25,gap = 5)
#     mkdashedstripes(alldata,r0 = 60,r1 = nrow-60,uu = [30,32,34],dur = 7,gap = 3)
#     mkdashedstripes(alldata,r0 = 90,r1 = nrow-90,uu = [50,52,54],dur = 5,gap = 2)
    if savetofile :
        np.savetxt("data03.txt",alldata,fmt = "%1d")
    return alldata
