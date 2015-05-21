
%module pysquid
%include "typemaps.i"
%include "carrays.i"
%array_class(double, dArray);
%{
#include "libsquid.h"
%}

// typemap to handle 4 element tuple output
%typemap(in,numinputs=0) double *TUPOUT_FOUR (double tmp[4]) {
   $1 = tmp;
}
%typemap(argout) double *TUPOUT_FOUR {
   PyObject *o = PyTuple_New(4);
   PyObject *o2, *o3; // temp values for appending to result
   // Put 4 output array values into Python Tuple
   int i;
   for (i=0; i<4; i++) {
      PyTuple_SetItem(o, i, PyFloat_FromDouble($1[i]));
   } 
   // Now append to $result if it already exists,
   // otherwise just return the tuple.
   if ((!$result)||($result == Py_None)) {
      $result = o;
   } else {
      if (!PyTuple_Check($result)) {
         PyObject *o2 = $result;
         $result = PyTuple_New(1);
         PyTuple_SetItem($result,0,o2);
      }
      o3 = PyTuple_New(1);
      PyTuple_SetItem(o3,0,o);
      o2 = $result;
      $result = PySequence_Concat(o2,o3);
      Py_DECREF(o2);
      Py_DECREF(o3);
   }
}

typedef unsigned long long squid_type;
// Note: typedef uint64_t squid_type doesn't work here!

#define DD2R 0.01745329252 // decimal deg to radians
#define THETAX 0.7297276562269663 // HSC proj polar dividing line (rad)
#define THETAXD 41.810314895778596 // HSC proj polar dividing line (deg)
#define KMAX 26 // maximum quad cube resolution

// Interpolation methods
// NEAREST = nearest neighbor
// BILINEAR = bilinear interpolation (2x2)
// CSPLINE = cubic b-spline interpolation (4x4)
// CCONVOL = cubic convolution interpolation (4x4)
enum {
  NEAREST, BILINEAR, CSPLINE, CCONVOL
};

// Map projections
// TSC = tangental spherical cube
// CSC = cobe spherical cube
// QSC = quadrilateralized spherical cube
// HSC = healpix spherical cube
enum {
  TSC, CSC, QSC, HSC
};

void sphdist(double lon1, double lat1, double lon2, double lat2, double *OUTPUT); // output=sdist
int face_range(int face, double x, double y, int *OUTPUT, double *OUTPUT, double *OUTPUT); // output=(newface,newx,newy)
int squid_validate(squid_type squid);
int squid_getres(squid_type squid);
int squid2xyfk(squid_type squid, squid_type *OUTPUT, squid_type *OUTPUT, int *OUTPUT, int *OUTPUT); // output=(x,y,face,k)
int xyfk2squid(squid_type x, squid_type y, int face, int k, squid_type *OUTPUT); // output=squid
int squid_getface(squid_type squid);
int squid_tside(int projection, int k, double cdelt, double *OUTPUT); // output=tside
int xyf2sph(int proj, double x, double y, int face, double *OUTPUT, double *OUTPUT); // output=(lon,lat)
int sph2xyf(int proj, double lon, double lat, double *OUTPUT, double *OUTPUT, int *OUTPUT); // output=(x,y,face)
int squid2sph(int proj, squid_type squid, double *OUTPUT, double *OUTPUT); // output=(lon,lat)
int sph2squid(int proj, double lon, double lat, int k, squid_type *OUTPUT); // output=squid
int squid_corners(int proj, squid_type squid, double *TUPOUT_FOUR, double *TUPOUT_FOUR); // output=(lon[4],lat[4])
int tile_xy2sph(int proj, squid_type squid, double x, double y, squid_type nside, double *OUTPUT, double *OUTPUT); //output=(lon,lat)
int tile_sph2xy(int proj, squid_type squid, double lon, double lat, squid_type nside, double *OUTPUT, double *OUTPUT);//output=(x,y)
int tile_nearest(int proj, squid_type squid, double lons, double lats, double *OUTPUT, double *OUTPUT); //output=(lon,lat)
void quadcube_getface(double lon, double lat, int *OUTPUT); // output=face
void hsc_getface(double lon, double lat, int *OUTPUT); // output=face

//--------------------------------------------------------------------
// cone_search() in python
//--------------------------------------------------------------------
%pythoncode %{

# Floating point version of python range function
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."
    if end == None:
        end = start + 0.0
        start = 0.0
    if inc == None:
        inc = 1.0
    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
    return L

#
# Cone search; lon,lat,rad are in radians
# kmax,kmin are max and min squid levels
# 
def cone_search(proj, lon, lat, rad, kmin, kmax):
   if (kmax > _pysquid.KMAX): # limit of squid resolution
      kmax = _pysquid.KMAX
   if (kmin > kmax): # make sure args are in right order
      tmp=kmax
      kmax=kmin
      kmin=tmp
   squidfull=[] # full tiles
   squidpart=[] # partial tiles
   # First check 6 face tiles
   (retval,squid0) = _pysquid.sph2squid(proj,lon,lat,0)
   for sp in range(8,14):
      # First check if all 4 corners are in search radius.
      # If so, it should be a fully contained tile.
      (retval,lonc,latc) = _pysquid.squid_corners(proj, sp)
      ccount=0 # corner count
      for i in range(0,4):
          cdist = _pysquid.sphdist(lonc[i],latc[i],lon,lat)
          if (cdist < rad):
             ccount=ccount+1
      if (ccount == 4): # all 4 points are in!
         if (kmin == 0):
            squidfull.append(sp)
         else:
            squidpart.append(sp)
         continue
      # Next check if search point is within tile
      if (sp == squid0):
         squidpart.append(sp)
         continue
      # Finally check if search region overlaps tile
      (retval,lon1,lat1) = _pysquid.tile_nearest(proj,sp,lon,lat)
      sdist = _pysquid.sphdist(lon1,lat1,lon,lat)
      if (sdist < rad):
         squidpart.append(sp)

   # Loop over resolution levels
   for k in range(0,kmax):
      # Get squid of lon,lat at this resolution level
      (retval,squid0) = _pysquid.sph2squid(proj,lon,lat,k+1)
      squidpart1=[] # partial squids at new resolution level
      # Now drop down a levela nd check subtiles of each
      # partial squid.
      for sp in squidpart:
         sp0=sp<<2 # shift over
         for i in range(0,4):
            sp1=sp0+i # get sub squid number
            # First check if all 4 corners are in search radius.
            # If so, it should be a fully contained tile.
            (retval,lonc,latc) = _pysquid.squid_corners(proj, sp1)
            ccount=0 # corner count
            for i in range(0,4):
               cdist = _pysquid.sphdist(lonc[i],latc[i],lon,lat)
               if (cdist < rad):
                  ccount=ccount+1
            if (ccount == 4): # all 4 points are in!
               if (kmin <= k+1):
                  squidfull.append(sp1)
               else:
                  squidpart1.append(sp1)
               continue
            # Next check if search point is within tile
            if (sp1 == squid0):
               squidpart1.append(sp1)
               continue
            # Finally check if search region overlaps tile
            (retval,lon1,lat1) = _pysquid.tile_nearest(proj,sp1,lon,lat)
            sdist = _pysquid.sphdist(lon1,lat1,lon,lat)
            if (sdist < rad):
               squidpart1.append(sp1)

      # Now update partial list
      squidpart=squidpart1
   
   return(squidfull,squidpart)
%}
