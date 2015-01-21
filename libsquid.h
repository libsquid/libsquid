//
// Library for performing Spherical-cube Quad-tree Unique ID (SQUID) computations.
// There are four possible projections: Tangental (TSC), COBE (CSC), Quad-Cube (QSC),
// and Healpix (HSC).
//
// -------------------------- LICENSE -----------------------------------
//
// This file is part of the LibSQUID software libraray.
//
// LibSQUID is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// LibSQUID is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright 2014 James Wren and Los Alamos National Laboratory
//

#ifndef LIBSQUID_H
#define LIBSQUID_H

#define LIBSQUID_VERSION "0.5.2"
#define LIBSQUID_MAJOR 0
#define LIBSQUID_MINOR 5
#define LIBSQUID_AGE 2
#define LIBSQUID_RELEASE ""

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <errno.h>

typedef uint64_t squid_type;

#ifndef DD2R // decimal degrees to radians
#define DD2R 0.01745329252
#endif

#ifndef PI
#define PI 3.141592653589793
#endif

#ifndef TWOPI
#define TWOPI 6.283185307179586
#endif

#ifndef HALFPI
#define HALFPI 1.5707963267948966
#endif

#ifndef THREEHALFPI
#define THREEHALFPI 4.7123889803846897
#endif

#ifndef QUARTPI
#define QUARTPI 0.7853981633974483
#endif

#ifndef THETAX // HSC lat dividing line, asin(2/3)
#define THETAX 0.7297276562269663
#define THETAXD 41.810314895778596
#endif

#ifndef LON_DIR // Longitude direction. 1 or -1
//#define LON_DIR 1 // longitude increases with x (common for geocentric coords)
#define LON_DIR -1 // longitude decreases with x (common for celestial coords)
#endif

#ifndef LON_POLE // Longitude of pole in rad, centerline of face1
//#define LON_POLE (-0.78539816340) // -45 deg
//#define LON_POLE 0.78539816340 // 45 deg
#define LON_POLE 0.0 // 0 deg
#endif

#ifndef KMAX // maximum quad cube resolution
#define KMAX 26
#endif

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

// max number of elements in squid array
#define ARR_MAX 1000

// max filename character length
#define FNAME_MAX 200

// set up utdate structure
struct ut_date {
  char *utdate;
  char *date;
  int year;
  int month;
  int day;
  double hour;
  double min;
  double sec;
  double mjd; // modified julian day
  time_t utime; // unix time
};

int utdate_populate(char *utdate, struct ut_date *datestr);
int utdate_future(struct ut_date *datestr, double duration, struct ut_date *futurestr);
void sph2rec(double lon, double lat, double *v);
void rec2sph(double *v, double *lon, double *lat);
void cross_prod(double *v1, double *v2, double *v3);
void sphdist(double lon1, double lat1, double lon2, double lat2, double *sdist);
void roty(double *v1, double ang, double *v2);
void rotz(double *v1, double ang, double *v2);
void arc_intercept(double r1, double d1, double r2, double d2, double r3, double d3,
		   double rc, double dc, double *r0, double *d0);
int interp_bilinear(double x, double y, double *img, long naxis1, long naxis2, double *outpix);
double bspline(double x);
double cubicon(double x);
int interp_bicubic(int method, double x, double y, double *img, long naxis1, long naxis2, double *outpix);
int interp_img(int method, double x, double y, double *img, long naxis1, long naxis2, double *outpix);
int qtree_id2xy(squid_type qtid, int k, squid_type *x, squid_type *y);
int qtree_xy2id(squid_type x, squid_type y, int k, squid_type *qtid);
int quadcube_arrange(int face, int *fnbr, int *frot);
int face_range(int face, double x, double y, int *newface, double *newx, double *newy);
int squid_validate(squid_type squid);
int squid_getres(squid_type squid);
int squid2xyfk(squid_type squid, squid_type *x, squid_type *y, int *face, int *k);
int xyfk2squid(squid_type x, squid_type y, int face, int k, squid_type *squid);
int squid_getface(squid_type squid);
int squid_tside(int projection, int k, double cdelt, double *tside);
int sph2LMN(double lon, double lat, int *face, double *L, double *M, double *N);
int LMNf2sph(double L, double M, double N, int face, double *lon, double *lat);
int xyf2sph(int proj, double x, double y, int face, double *lon, double *lat);
int sph2xyf(int proj, double lon, double lat, double *x, double *y, int *face);
int tsc_xyf2sph(double x, double y, int face, double *lon, double *lat);
int tsc_sph2xyf(double lon, double lat, double *x, double *y, int *face);
double csc_fwarp(double x, double y);
double csc_rwarp(double x, double y);
int csc_xyf2sph(double x, double y, int face, double *lon, double *lat);
int csc_sph2xyf(double lon, double lat, double *x, double *y, int *face);
int qsc_xyf2sph(double x, double y, int face, double *lon, double *lat);
int qsc_sph2xyf(double lon, double lat, double *x, double *y, int *face);
int hsc_xyf2sph(double x, double y, int face, double *lon, double *lat);
int hsc_sph2xyf(double lon, double lat, double *x, double *y, int *face);
int squid2sph(int proj, squid_type squid, double *lon, double *lat);
int sph2squid(int proj, double lon, double lat, int k, squid_type *squid);
int squid_corners(int proj, squid_type squid, double *lonrr, double *darr);
int tile_xy2sph(int proj, squid_type squid, double x, double y, squid_type tside, double *lon, double *lat);
int tile_sph2xy(int proj, squid_type squid, double lon, double lat, squid_type tside, double *x, double *y);
int tile_nearest(int proj, squid_type squid, double lon, double lats, double *lonn, double *latn);
int cone_search(int proj, double lon, double lat, double srad, int kmin, int kmax,
		long *nfull, squid_type **full_tiles, long *npart, squid_type **part_tiles);
void quadcube_getface(double lon, double lat, int *face);
void hsc_getface(double lon, double lat, int *face);

#ifdef __cplusplus
}
#endif

#endif //LIBSQUID_H

