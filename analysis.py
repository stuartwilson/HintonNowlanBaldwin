import numpy as np
import pylab as pl
import h5py
import sys

if (len(sys.argv)<2):
    fname = 'logs/out.h5'
else:
    fname = sys.argv[1]

F = h5py.File(fname)
p0 = F['p0'][:]
p1 = F['p1'][:]
p2 = F['p2'][:]
pCorr = F['pCorr'][:]
pIncorr = 1-pCorr-p2

F.close()

'''
F = pl.figure()
f = F.add_subplot(211)
f.plot(p0)
f.plot(p1)
f.plot(p2)
f.legend(['0','1','?'])

f = F.add_subplot(212)
f.plot(pIncorr)
f.plot(pCorr)
f.plot(p2)
f.legend(['Incorr','Corr','?'])
'''



## TRAJECTORY PLOT

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # appropriate import to draw 3d polygons
from matplotlib import style

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

F = plt.figure(figsize=(9,9))

gr = 0.7
len = 1.1
greeny = (0.5,0.6,0.5,0.5)

#####
f = F.add_subplot(221,projection='3d')
f.plot3D([1,0],[0,0],[0,1],'--',color=(0,0,0),linewidth=1)
f.plot3D([0,0],[0,1],[1,0],'-',color=(0,0,1),linewidth=1)
f.plot3D([0,0],[0,0],'-',color=greeny,linewidth=8)
f.plot3D([0,0],[0,0],[0,len],color=(gr,gr,gr),linewidth=1)
f.plot3D([0,0],[len,0],[0,0],color=(gr,gr,gr),linewidth=1)
f.plot3D([0,len],[0,0],[0,0],color=(gr,gr,gr),linewidth=1)

f.plot3D(p0,p2,p1,'-',color=(1,0,0),linewidth=2)
f.scatter(p0[-1],p2[-1],p1[-1],color=(0,0,0))
len2 = len+0.1
f.text(len2,0,0,r'$p(0)$')
f.text(0,len2,0,r'$p(?)$')
f.text(0,0,len2,r'$p(1)$')
#x-2y+z=6
X1=np.array([0, 0, 1])
Y1=np.array([0, 1, 0])
Z1=np.array([1, 0, 0])
verts = [list(zip(X1, Y1, Z1))]
srf = Poly3DCollection(verts, alpha=greeny[3], facecolor=(greeny[0],greeny[1],greeny[2]))
plt.gca().add_collection3d(srf)
axisEqual3D(f)
f.axes.get_xaxis().set_visible(False)
f.axes.get_yaxis().set_visible(False)
f.axes.get_zaxis().set_visible(False)
f.set_xticks([])
f.set_yticks([])
f.set_zticks([])
f.w_xaxis.set_pane_color((1.,1.,1.,1.))
f.w_yaxis.set_pane_color((1.,1.,1.,1.))
f.w_zaxis.set_pane_color((1.,1.,1.,1.))
#f.legend([r'$p(1)=1-p(0)$',r'$p(1)=1-p(?)$','valid'],frameon=False)

f.set_axisbelow(True)
f.grid('off')
f.view_init(elev=56, azim=18)
f.set_xlabel('X')
f.set_ylabel('Y')
f.set_zlabel('Z')
f.set_axis_off()

#####
f = F.add_subplot(222,projection='3d')
f.plot3D([1,0],[0,0],[0,1],'--',color=(0,0,0),linewidth=1)
f.plot3D([0,0],[0,1],[1,0],'-',color=(0,0,1),linewidth=1)
f.plot3D([0,0],[0,0],'-',color=greeny,linewidth=8)
f.plot3D([0,0],[0,0],[0,len],color=(gr,gr,gr),linewidth=1)
f.plot3D([0,0],[len,0],[0,0],color=(gr,gr,gr),linewidth=1)
f.plot3D([0,len],[0,0],[0,0],color=(gr,gr,gr),linewidth=1)

f.plot3D(pIncorr,p2,pCorr,'-',color=(1,0,0),linewidth=2)
f.scatter(pIncorr[-1],p2[-1],pCorr[-1],color=(0,0,0))
len2 = len+0.1
f.text(len2,0,0,r'$p(N)$')
f.text(0,len2,0,r'$p(?)$')
f.text(0,0,len2,r'$p(Y)$')
#x-2y+z=6
X1=np.array([0, 0, 1])
Y1=np.array([0, 1, 0])
Z1=np.array([1, 0, 0])
verts = [list(zip(X1, Y1, Z1))]
srf = Poly3DCollection(verts, alpha=greeny[3], facecolor=(greeny[0],greeny[1],greeny[2]))
plt.gca().add_collection3d(srf)
axisEqual3D(f)
f.axes.get_xaxis().set_visible(False)
f.axes.get_yaxis().set_visible(False)
f.axes.get_zaxis().set_visible(False)
f.set_xticks([])
f.set_yticks([])
f.set_zticks([])
f.w_xaxis.set_pane_color((1.,1.,1.,1.))
f.w_yaxis.set_pane_color((1.,1.,1.,1.))
f.w_zaxis.set_pane_color((1.,1.,1.,1.))
#f.legend([r'$p(1)=1-p(0)$',r'$p(1)=1-p(?)$','valid'],frameon=False)

f.set_axisbelow(True)
f.grid('off')
f.view_init(elev=56, azim=18)
f.set_xlabel('X')
f.set_ylabel('Y')
f.set_zlabel('Z')
f.set_axis_off()


f = F.add_subplot(223)
f.plot(p0)
f.plot(p1)
f.plot(p2)
f.legend(['0','1','?'])

f = F.add_subplot(224)
f.plot(pIncorr)
f.plot(pCorr)
f.plot(p2)
f.legend(['N','Y','?'])


F.tight_layout()
plt.show()


pl.show()
