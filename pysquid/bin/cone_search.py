#!/usr/bin/python -tt

# -------------------------- LICENSE -----------------------------------
#
# This file is part of the LibSQUID software libraray.
#
# LibSQUID is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LibSQUID is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014 James Wren and Los Alamos National Laboratory
#

import os
import sys
import math as m
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../")
import pysquid as ps

# constants
dd2r=ps.DD2R # decimal deg to radians

# get args
if (len(sys.argv) < 7):
    print "Find squid tiles within search region."
    print "calling sequence: cone_test.py projection lon lat rad kmin kmax"
    print "projections: 0=TSC, 1=CSC, 2=QSC, 3=HSC"
    print "kmin,kmax are min and max squid levels for search"
    print "angles in degrees"
    sys.exit()
proj=int(sys.argv[1])
ra=float(sys.argv[2])*dd2r # lon
dec=float(sys.argv[3])*dd2r # lat
rad=float(sys.argv[4])*dd2r
kmin=int(sys.argv[5])
kmax=int(sys.argv[6])

# check projection
if (proj == 0):
    print "using TSC projection"
elif (proj == 1):
    print "using CSC projection"
elif (proj == 2):
    print "using QSC projection"
else:
    print "using HSC projection"
    proj=3

# get appropriate squid tiles for search region
(fullsquids,partsquids)=ps.cone_search(proj,ra,dec,rad,kmin,kmax)
print "List of full tiles in search region:"
print fullsquids
print "List of partial tiles in search region:"
print partsquids

# function to convert face x,y to world x,y
def GetWorldXY(x,y,f):
    if (int(f)==0):
        xw=x+2.0
        yw=y+2.0
    elif (int(f)==5):
        xw=x+2.0
        yw=y
    elif (int(f)==4):
        xw=x+3.0
        yw=y
    else:
        xw=x+(3.0-f)
        yw=y+1.0
    return(xw,yw)

# Spherical to Cartesian conversion
# ra,dec in radians
def sph2rec(ra,dec):
    phi=ra
    theta=m.pi/2.0-dec
    x=m.sin(theta)*m.cos(phi)
    y=m.sin(theta)*m.sin(phi)
    z=m.cos(theta)
    return (x,y,z)

# Cartesian to Spherical converstion
# ra,dec output in radians
def rec2sph(x,y,z):
    r=m.sqrt(x**2+y**2+z**2)
    theta=m.acos(z/r)
    phi=m.atan2(y,x)
    ra=phi
    dec=m.pi/2.0-theta
    return (ra,dec)

# Rotate about y axis (phi=90,theta=90) ang in rad
def roty(x,y,z,ang):
    arr1=[m.cos(ang),0.0,-1.0*m.sin(ang)]
    arr2=[0.0,1.0,0.0]
    arr3=[m.sin(ang),0.0,m.cos(ang)]
    Ry=np.matrix([arr1,arr2,arr3])
    cart=np.matrix([[x,y,z]])
    cart2=Ry*cart.transpose()
    return cart2[0,0],cart2[1,0],cart2[2,0]

# Rotate about z axis (phi=0,theta=0) ang in rad
def rotz(x,y,z,ang):
    arr1=[m.cos(ang),-1.0*m.sin(ang),0.0]
    arr2=[m.sin(ang),m.cos(ang),0.0]
    arr3=[0.0,0.0,1.0]
    Rz=np.matrix([arr1,arr2,arr3])
    cart=np.matrix([[x,y,z]])
    cart2=Rz*cart.transpose()
    return cart2[0,0],cart2[1,0],cart2[2,0]

#
# Plot search region
#
(retval,xc,yc,fc)=ps.sph2xyf(proj,ra,dec)
(xcw,ycw)=GetWorldXY(xc,yc,fc)
hx=[]
hy=[]
dang=m.pi/2.0-dec
for rax in ps.frange(0,2.0*m.pi,.001):
    decx=m.pi/2.0-rad
    (x0,y0,z0)=sph2rec(rax,decx)
    (x1,y1,z1)=roty(x0,y0,z0,dang)
    (x2,y2,z2)=rotz(x1,y1,z1,ra-m.pi)
    (ra1,dec1)=rec2sph(x2,y2,z2)
    (retval,hx0,hy0,hf)=ps.sph2xyf(proj,ra1,dec1)
    (hxw,hyw)=GetWorldXY(hx0,hy0,hf)
    hx.append(hxw)
    hy.append(hyw)
plt.plot([xcw],[ycw],'ro')
plt.plot(hx,hy,'r,')

#
# Loop over partial tiles and plot them
#
minpx=1e6
maxpx=-1e6
minpy=1e6
maxpy=-1e6
for squid in partsquids:
    px=[]
    py=[]
    (retval,rac,decc)=ps.squid_corners(proj,squid)
    for i in [0,1,3,2,0]:
        (retval,px0,py0,pf)=ps.sph2xyf(proj,rac[i],decc[i])
        (pxw,pyw)=GetWorldXY(px0,py0,pf)
        if (pxw < minpx):
            minpx=pxw
        if (pxw > maxpx):
            maxpx=pxw
        if (pyw < minpy):
            minpy=pyw
        if (pyw > maxpy):
            maxpy=pyw
        px.append(pxw)
        py.append(pyw)
    plt.plot(px,py,'b')
    stxt=" {:d}".format(squid)
    plt.text(px[0],py[0],stxt)
    

#
# Loop over full tiles and plot them
#
minfx=1e6
maxfx=-1e6
minfy=1e6
maxfy=-1e6
for squid in fullsquids:
    fx=[]
    fy=[]
    (retval,rac,decc)=ps.squid_corners(proj,squid)
    for i in [0,1,3,2,0]:
        (retval,fx0,fy0,ff)=ps.sph2xyf(proj,rac[i],decc[i])
        (fxw,fyw)=GetWorldXY(fx0,fy0,ff)
        if (fxw < minfx):
            minfx=fxw
        if (fxw > maxfx):
            maxfx=fxw
        if (fyw < minfy):
            minfy=fyw
        if (fyw > maxfy):
            maxfy=fyw
        fx.append(fxw)
        fy.append(fyw)
    plt.plot(fx,fy,'g')
    stxt=" {:d}".format(squid)
    plt.text(fx[0],fy[0],stxt,color='g')

minx=min([min(hx),minpx,minfx])
maxx=max([max(hx),maxpx,maxfx])
miny=min([min(hy),minpy,minfy])
maxy=max([max(hy),maxpy,maxfy])
plt.axis([minx-.01,maxx+.01,miny-.01,maxy+.01])
plt.show()
