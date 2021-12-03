import numpy as np
import h5py
import sys
import os

cutoffIncrement = 10
P0 = np.zeros([100,1000])
P1 = np.zeros([100,1000])
P2 = np.zeros([100,1000])

a = os.listdir('logs')
for i in a:
    b=os.listdir('logs/'+i)
    cutoff = int(i.split('ff')[1])      # split at ff in e.g., cutoff20
    cutoffInd = int(cutoff/cutoffIncrement)-1        #
    for j in b:
        simInd = int(j.split('m')[1])-1 # split at m in e.g., sim57
        fname = 'logs/'+i+'/'+j+'/out.h5'

        Fl = h5py.File(fname)
        p0 = Fl['p0'][:]
        p1 = Fl['p1'][:]
        p2 = Fl['p2'][:]

        Fl.close()

        P0[cutoffInd,simInd] = p0[-1]
        P1[cutoffInd,simInd] = p1[-1]
        P2[cutoffInd,simInd] = p2[-1]

np.savez('data.npz',P0=P0,P1=P1,P2=P2)
