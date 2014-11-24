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

dd2r=ps.DD2R

# get args (in deg)
if (len(sys.argv) < 4):
    print "Determine the closest point on squid tile to given coords."
    print "calling sequence: ./tile_nearest.py projection squid lon lat"
    print "projections; 0=TSC, 1=CSC, 2=QSC, 3=HSC"
    print "lon,lat in decimal degrees"
    sys.exit()
proj=int(sys.argv[1])
squid=long(sys.argv[2])
ra=float(sys.argv[3])*dd2r # lon
dec=float(sys.argv[4])*dd2r # lat

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

# Get nearest point on squid tile border to ra,dec
(rval,ra0,dec0)=ps.tile_nearest(proj,squid,ra,dec)
print "nearest lon={:f} lat={:f} (deg)".format(ra0/dd2r,dec0/dd2r)
rad=ps.sphdist(ra0,dec0,ra,dec)
print "sphdist={:f} (deg)".format(rad/dd2r)

# Get params for this squid
sk=ps.squid_getres(squid) # squid resolution parameter
(rval,sr,sd)=ps.squid2sph(proj,squid) # ra,dec of squid center
(rval,sx,sy,sf)=ps.sph2xyf(proj,sr,sd) # get x,y,face for squid center

# cone array
(rval,hxc,hyc,hfc)=ps.sph2xyf(proj,ra,dec) # x,y,face for search point
(rval,xp0,yp0,pf)=ps.sph2xyf(proj,ra0,dec0) # x,y,face for nearest point

# get array for search radius for plotting
hx=[]
hy=[]
dang=m.pi/2.0-dec
for rax in ps.frange(0,2.0*m.pi,.001):
    decx=m.pi/2.0-rad
    (x0,y0,z0)=sph2rec(rax,decx)
    (x1,y1,z1)=roty(x0,y0,z0,dang)
    (x2,y2,z2)=rotz(x1,y1,z1,ra-m.pi)
    (ra1,dec1)=rec2sph(x2,y2,z2)
    (rval,hx0,hy0,hf)=ps.sph2xyf(proj,ra1,dec1)
    if (hf == sf):
        hx.append(hx0)
        hy.append(hy0)

# plot gridlines
nlines=int(m.pow(2,sk))
plt.axis([-.1,1.1,-.1,1.1])
for i in range(0,nlines+1):
    r=1.0*i/nlines
    plt.plot([r,r],[0,1],'b')
    plt.plot([0,1],[r,r],'b')

if (hfc == sf):
    plt.plot([hxc],[hyc],'ro')
plt.plot(hx,hy,'r,')
plt.plot(sx,sy,'go')
if (pf == sf):
    plt.plot([xp0],[yp0],'gx')
else:
    plt.plot([xp0],[yp0],'rx')
plt.show()
