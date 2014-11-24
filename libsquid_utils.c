//
// Utility functions for libsquid library.
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

#include <libsquid.h>

// For some reason regular "round" causes problems.
#ifndef round
double round(double round);
/* Based on http://forums.belution.com/en/cpp/000/050/13.shtml */
double round(double value)
{
    if (value < 0)
        return -(floor(-value + 0.5));
    else
        return   floor( value + 0.5);
}
#endif

// Populate utdate struct given utdate string
int utdate_populate(char *utdate, struct ut_date *datestr) {
  char *utime, uyr[3], umn[3], udy[3], *uhr, *umin;
  char *tokptr;
  struct tm uttm;

  datestr->utdate=strdup(utdate);
  datestr->date=strtok_r(utdate,"T",&utime);
  if ((datestr->date == NULL)||(strlen(datestr->date) != 6)) {
    fprintf(stderr,"improper date format in utdate_populate\n");
    return(-1);
  }
  uyr[0]=datestr->date[0]; uyr[1]=datestr->date[1]; uyr[2]='\0';
  umn[0]=datestr->date[2]; umn[1]=datestr->date[3]; umn[2]='\0';
  udy[0]=datestr->date[4]; udy[1]=datestr->date[5]; udy[2]='\0';
  datestr->year=atoi(uyr);
  datestr->month=atoi(umn);
  datestr->day=atoi(udy);
  uhr=strtok_r(utime,":",&tokptr);
  if (uhr != NULL) {
    datestr->hour=atof(uhr);
    utime=tokptr;
    umin=strtok_r(utime,":",&tokptr);
    if (umin != NULL) {
      datestr->min=atof(umin);
      if (tokptr != NULL) datestr->sec=atof(tokptr);
    }
  }

  uttm.tm_sec=(int)datestr->sec;
  uttm.tm_min=(int)datestr->min;
  uttm.tm_hour=(int)datestr->hour;
  uttm.tm_mday=datestr->day;
  uttm.tm_mon=datestr->month-1;
  uttm.tm_year=datestr->year+100;
  uttm.tm_isdst=0;
  if ((datestr->utime=timegm(&uttm))==-1) {
    fprintf(stderr,"timegm failed in utdate_populate\n");
    return(-1);
  }

  // mjd of 1/1/1970 0 UT was 40587.00
  datestr->mjd=(double)datestr->utime/86400.0+40587.00;

  return(0);
}

// Populate utdate struct given utdate string for future time.
// Duration is in decimal hours
int utdate_future(struct ut_date *datestr, double duration, struct ut_date *futurestr) {
  struct tm *uttm;
  char utdate[FNAME_MAX];
  char date[FNAME_MAX];
  
  futurestr->utime=datestr->utime+(time_t)(duration*3600.0);
  if ((uttm=gmtime(&futurestr->utime))==NULL) {
    fprintf(stderr,"gmtime failed in utdate_future\n");
    return(-1);
  }
  futurestr->sec=uttm->tm_sec;
  futurestr->min=uttm->tm_min;
  futurestr->hour=uttm->tm_hour;
  futurestr->day=uttm->tm_mday;
  futurestr->month=uttm->tm_mon+1;
  futurestr->year=uttm->tm_year-100;
  futurestr->mjd=datestr->mjd+(duration/24.0);

  snprintf(date,FNAME_MAX-1,"%02d%02d%02d",futurestr->year, futurestr->month, futurestr->day);
  futurestr->date=strdup(date);
  snprintf(utdate,FNAME_MAX-1,"%02d%02d%02dT%02d:%02d:%05.2f",
		futurestr->year, futurestr->month, futurestr->day,
		(int)futurestr->hour,(int)futurestr->min,futurestr->sec);
  futurestr->utdate=strdup(utdate);
  

  return(0);
}

// Spherical to Cartesian conversion
// lon,lat in radians
// v = 3 element output vector [x,y,z]
void sph2rec(double lon, double lat, double *v) {
  double theta;
  theta=HALFPI-lat;
  v[0]=sin(theta)*cos(lon);
  v[1]=sin(theta)*sin(lon);
  v[2]=cos(theta);
}
  
// Cartesian to Spherical conversion
// lon,lat in radians
void rec2sph(double *v, double *lon, double *lat) {
  double r, theta;
  r=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  theta=acos(v[2]/r);
  *lon=atan2(v[1],v[0]);
  *lat=HALFPI-theta;
}

// Cross product
// v1 and v2 are 3 element input vectors
// v3 is 3 element output vector
void cross_prod(double *v1, double *v2, double *v3) {
  v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
  v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
  v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

// Spherical distance between two points
// lon,lat in radians
void sphdist(double lon1, double lat1, double lon2, double lat2, double *sdist) {
  double v1[3], v2[3];
  double d,xc,yc,zc,sn;
  sph2rec(lon1,lat1,v1);
  sph2rec(lon2,lat2,v2);
  d=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  xc=v1[1]*v2[2]-v1[2]*v2[1];
  yc=v1[2]*v2[0]-v1[0]*v2[2];
  zc=v1[0]*v2[1]-v1[1]*v2[0];
  sn=sqrt(xc*xc+yc*yc+zc*zc);
  *sdist=atan2(sn,d);
}

// rotate about y axis (phi=90, theta=90)
// v1 = 3 element cartesian input vector
// ang = rotation angle in rad
// v2 = 3 element cartesian output vector
void roty(double *v1, double ang, double *v2) {
  double m1x, m1y, m1z, m2x, m2y, m2z, m3x, m3y, m3z; // rotation matrix
  m1x=cos(ang); m1y=0.0; m1z=-1.0*sin(ang);
  m2x=0.0; m2y=1.0; m2z=0.0;
  m3x=sin(ang); m3y=0.0; m3z=cos(ang);
  v2[0]=m1x*v1[0]+m1y*v1[1]+m1z*v1[2];
  v2[1]=m2x*v1[0]+m2y*v1[1]+m2z*v1[2];
  v2[2]=m3x*v1[0]+m3y*v1[1]+m3z*v1[2];
}

// rotate about z axis (phi=0, theta=0)
// v1 = 3 element cartesian input vector
// ang = rotation angle in rad
// v2 = 3 element cartesian output vector
void rotz(double *v1, double ang, double *v2) {
  double m1x, m1y, m1z, m2x, m2y, m2z, m3x, m3y, m3z; // rotation matrix
  m1x=cos(ang); m1y=-1.0*sin(ang); m1z=0.0;
  m2x=sin(ang); m2y=cos(ang); m2z=0.0;
  m3x=0.0; m3y=0.0; m3z=1.0;
  v2[0]=m1x*v1[0]+m1y*v1[1]+m1z*v1[2];
  v2[1]=m2x*v1[0]+m2y*v1[1]+m2z*v1[2];
  v2[2]=m3x*v1[0]+m3y*v1[1]+m3z*v1[2];
}

// Given 3 points on a sphere defininig an arc and 4th point (rc,dc),
// find nearest point on the arc to the 4th point.  r2,d2 should be
// the point between r1,d1 and r3,d3.  Arc points must be
// coplanar or consider this an estimate.  Angles in radians.
// If the intercept is outside the range of the arc, the
// nearest end is chosen.  r0,d0 are the coords of the output point
// Note: here r? is lon and d? is lat.
void arc_intercept(double r1, double d1, double r2, double d2, double r3, double d3,
		   double rc, double dc, double *r0, double *d0) {
  double cart1[3], cart2[3], cart3[3], carts[3]; // cartesian coords for r1,d1 r2,d2 r3,d3 and rs,ds
  double v1[3], v2[3]; // vectors r1,d1->r2,d2 and r2,d2->r3,d3
  double F[3]; // cross product of v1xv2
  double lonf,latf; // sph coords of F
  double ang1, ang2, ang3; // rotation angles in rad
  double vr1[3], vr2[3], vr3[3]; // temporary rotate vectors
  double rlonm0, rlatm0, rlonm,rlatm; // lon,lat of midpoint after F rotated to pole
  double rlon1, rlat1, rlon3, rlat3, rlons, rlats; // end points and search point rotated
  double minrlon, maxrlon; // min([rlon1,rlon2]) and max([rlon1,rlon2])
  double vn[3]; // nearest point in cartesian
  double lon0,lat0; // nearest point in r,d
  double sdist1, sdist2; // sph dist of endpoints to rc,dc

  // Convert all three arc point to cartesian and get vectors
  sph2rec(r1,d1,cart1);
  sph2rec(r2,d2,cart2);
  sph2rec(r3,d3,cart3);
  v1[0]=cart2[0]-cart1[0];
  v1[1]=cart2[1]-cart1[1];
  v1[2]=cart2[2]-cart1[2];
  v2[0]=cart3[0]-cart2[0];
  v2[1]=cart3[1]-cart2[1];
  v2[2]=cart3[2]-cart2[2];
  // Get cross product
  cross_prod(v1,v2,F);
  rec2sph(F,&lonf,&latf);
  ang1=2*PI-lonf;
  ang2=HALFPI-latf;

  // Get midpoint coords after rotating F to pole
  rotz(cart2,ang1,vr1);
  roty(vr1,ang2,vr2);
  rec2sph(vr2,&rlonm0,&rlatm0);
  // rotate one more time to put midpoint at 180
  ang3=PI-rlonm0;
  rotz(vr2,ang3,vr3);
  rec2sph(vr3,&rlonm,&rlatm);
  rlonm=fmod(rlonm+2*PI,2*PI);

  // Rotate search point in same way
  sph2rec(rc,dc,carts);
  rotz(carts,ang1,vr1);
  roty(vr1,ang2,vr2);
  rotz(vr2,ang3,vr3);
  rec2sph(vr3,&rlons,&rlats);
  rlons=fmod(rlons+2*PI,2*PI);
  // Rotate end point #1 in same way
  rotz(cart1,ang1,vr1);
  roty(vr1,ang2,vr2);
  rotz(vr2,ang3,vr3);
  rec2sph(vr3,&rlon1,&rlat1);
  rlon1=fmod(rlon1+2*PI,2*PI);
  // Rotate end point #2 in same way
  rotz(cart3,ang1,vr1);
  roty(vr1,ang2,vr2);
  rotz(vr2,ang3,vr3);
  rec2sph(vr3,&rlon3,&rlat3);
  rlon3=fmod(rlon3+2*PI,2*PI);

  // If rlons inside range of endpoints we're good!
  minrlon = (rlon1 < rlon3) ? rlon1 : rlon3;
  maxrlon = (rlon1 > rlon3) ? rlon1 : rlon3;
  if ((rlons < maxrlon)&&(rlons > minrlon)) {
    // Get new point and rotate back
    sph2rec(rlons,rlatm,vn);
    rotz(vn,-1.0*ang3,vr1);
    roty(vr1,-1.0*ang2,vr2);
    rotz(vr2,-1.0*ang1,vr3);
    rec2sph(vr3,&lon0,&lat0);
    lon0=fmod(lon0+2*PI,2*PI);
    *r0=lon0;
    *d0=lat0;
    return;
  }

  // If we got here, it must be one of the end points.
  // Get closest endpoint and return that.
  sphdist(r1,d1,rc,dc,&sdist1);
  sphdist(r3,d3,rc,dc,&sdist2);
  if (sdist1 < sdist2) {
    *r0=r1;
    *d0=d1;
  } else {
    *r0=r3;
    *d0=d3;
  }
  return;
}

// Bilinear interpolation
// x,y are the actual coordinates to be interpolated
// pixel at img[0,0] has center x,y of [0,0]
// naxis1,naxis2 is the width of your image x,y such that the last
//   pixel index is [naxis1-1,naxis2-1].
// outpix is the interpolated pixel value
int interp_bilinear(double x, double y, double *img, long naxis1, long naxis2, double *outpix) {
  long x0,y0,x1,y1; // pixel coords surrounding x,y
  double xd,yd; // distance from x,y to x0,y0
  double pix00,pix01,pix10,pix11; // pixel values surrounding x,y
  long pix; // image pixel location
  // Check for valid pixel value
  if ((x < 0)||(y < 0)||(x >= naxis1)||(y >= naxis2)) {
    fprintf(stderr,"pixel value out of bounds for bilinear interpolation\n");
    fprintf(stderr,"x=%.5f y=%.5f\n",x,y);
    *outpix=0.0;
    return(-1);
  }
  // Get adjacent pix coords
  x0=(long)floor(x);
  y0=(long)floor(y);
  x1=x0+1;
  y1=y0+1;
  xd=x-x0;
  yd=y-y0;
  // Handle edge pixels
  if (x0 < 0) x0=x1;
  if (x1 >= naxis1) x1=x0;
  if (y0 < 0) y0=y1;
  if (y1 >= naxis2) y1=y0;
  // Get adjacent pixel values
  pix=(long)(x0+(naxis1*y0));
  pix00=img[pix];
  pix=(long)(x0+(naxis1*y1));
  pix01=img[pix];
  pix=(long)(x1+(naxis1*y0));
  pix10=img[pix];
  pix=(long)(x1+(naxis1*y1));
  pix11=img[pix];
  // Calculate output pixel
  *outpix=((1-xd)*(((1-yd)*pix00)+(yd*pix01))) + (xd*(((1-yd)*pix10)+(yd*pix11)));
  return(0);
}

// B-spline interpolation method
double bspline(double x) {
  double xa, bval;
  xa=fabs(x);
  if ((xa >= 0)&&(xa <= 1)) {
    bval=(2.0/3.0)+(pow(xa,3.0)/2.0)-pow(x,2.0);
  } else {
    bval=pow(2.0-xa,3.0)/6.0;
  }
  return(bval);
}
// Cubic convolution interpolation method
double cubicon(double x) {
  double xa, bval;
  double ca=-0.5; // weighting factor
  xa=fabs(x);
  if ((xa >= 0)&&(xa <= 1)) {
    bval=((ca+2.0)*pow(xa,3.0))-((ca+3.0)*pow(xa,2.0))+1.0;
  } else {
    bval=(ca*pow(xa,3.0))-(5.0*ca*pow(xa,2.0))+(8.0*ca*xa)-(4.0*ca);
  }
  return(bval);
}
// Bi-cubic interpolation
// x,y are the actual coordinates to be interpolated
// pixel at img[0,0] has center x,y of [0,0].
// naxis1,naxis2 is the width of your image x,y such that the last
//   pixel index is [naxis1-1,naxis2-1].
// outpix is the interpolated pixel value.
// method is either CSPLINE (bspline) or CCONVOL (cubic convolution, default)
int interp_bicubic(int method, double x, double y, double *img, long naxis1, long naxis2, double *outpix) {
  double pixval; // temporary value for outpix
  double bvalx,bvaly; // temporary bspline values
  long x0,y0; // coords of nearest pix center below x,y
  double xd,yd; // distance from x,y to x0,y0
  long m,n; // loop index values
  long pix; // image pixel location
  x0=(long)floor(x);
  y0=(long)floor(y);
  xd=x-x0;
  yd=y-y0;
  pixval=0.0;
  // Check for edge conditions, if too close revert to bilinear
  if ((x < 2)||(x > (naxis1-3))||(y < 2)||(y > (naxis2-3))) {
    if (interp_bilinear(x,y,img,naxis1,naxis2,outpix) < 0) {
      fprintf(stderr,"bilinear interpolation failed in interp_bicubic\n");
      return(-1);
    }
    return(0);
  }
  // Do interpolation
  for (m=-1; m<3; m++) {
    for (n=-1; n<3; n++) {
      if (method == CSPLINE) {
	bvalx=bspline(m-xd);
	bvaly=bspline(-1.0*(n-yd));
      } else {
	bvalx=cubicon(m-xd);
	bvaly=cubicon(-1.0*(n-yd));
      }
      pix=(long)((x0+m)+(naxis1*(y0+n)));
      pixval=pixval+(img[pix]*bvalx*bvaly);
    }
  }
  *outpix=pixval;
  return(0);
}

// General interpolation call (uses one of the methods above)
int interp_img(int method, double x, double y, double *img, long naxis1, long naxis2, double *outpix) {
  double xr,yr;
  long pix;
  if (method == NEAREST) {
    xr=round(x);
    yr=round(y);
    if ((xr < 0)||(xr >= naxis1)||(yr < 0)||(yr >= naxis2)) {
      // pixel values out of range
      fprintf(stderr,"pixel value out of bounds in interp_img");
      *outpix=0.0;
      return(-1);
    }
    pix=(long)(xr+(naxis1*yr));
    *outpix=img[pix];
    return(0);
  } else if (method != BILINEAR) { // bicubic interpolation, either CSPLINE or CCONVOL
    // default here with be CCONVOL
    if (interp_bicubic(method,x,y,img,naxis1,naxis2,outpix) < 0) {
      fprintf(stderr,"interp_bicubic failed in interp_image\n");
      return(-1);
    }
    return(0);
  } else { // default BILINEAR
    if (interp_bilinear(x,y,img,naxis1,naxis2,outpix) < 0) {
      fprintf(stderr,"interp_bilinear failed in interp_image\n");
      return(-1);
    }
    return(0);
  }
  return(0); // redundant, but may be needed in future
}
    

// quad-tree ID to x,y
int qtree_id2xy(squid_type qtid, int k, squid_type *x, squid_type *y) {
  int bitlen, nbit;
  squid_type ix,iy,tmpbit;

  ix=0;
  iy=0;
  bitlen = 2*k;
  for (nbit=0; nbit < bitlen; nbit++) {
    tmpbit = (qtid>>nbit)&1;
    if (nbit%2) {
      iy=iy+(tmpbit<<((nbit-1)/2));
    } else {
      ix=ix+(tmpbit<<(nbit/2));
    }
  }
  *x=ix;
  *y=iy;
  return(0);
}

// x,y to quad-tree id
int qtree_xy2id(squid_type x, squid_type y, int k, squid_type *qtid) {
  int nbit;
  squid_type nside;
  squid_type xid,yid,tmpbit;

  nside=(squid_type)pow(2,k);
  if ((x > nside)|| (y > nside)) {
    // x,y out of range
    fprintf(stderr,"x or y out of range in qtree_xy2id\n");
    return(-1);
  }
  xid=0;
  yid=0;
  for (nbit=0; nbit < k; nbit++) {
    tmpbit=(x>>nbit)&1;
    xid=xid+(tmpbit<<(nbit*2));
    tmpbit=(y>>nbit)&1;
    yid=yid+(tmpbit<<(nbit*2+1));
  }
  *qtid=xid|yid;
  return(0);
}

// Quadcube face arrangement
int quadcube_arrange(int face, int *fnbr, int *frot) {
  // fnbr = [up,down,left,right,back] = face neighbors
  // lat increases upward, lon increases to the right (if LON_DIR=1)
  // frot=rotation array, ccw, angle in radians *2/pi
  // rot=0 >> no rotation x'=x, y'=y
  // rot=1 >> 90 deg ccw, x'=-y, y'=x
  // rot=2 >> 180 deg ccw, x'=-x, y'=-y
  // rot=3 >> 270 deg ccw, x'=y, y'=-x

  switch (face) {
  case 0:
    //fnbr=[3,1,2,4,5];
    fnbr[0]=3; fnbr[1]=1; fnbr[2]=2; fnbr[3]=4; fnbr[4]=5;
    //frot=[2,0,3,1,2];
    frot[0]=2; frot[1]=0; frot[2]=3; frot[3]=1; frot[4]=2;
    break;
  case 1:
    // the prototype orientation
    //fnbr=[0,5,2,4,3];
    fnbr[0]=0; fnbr[1]=5; fnbr[2]=2; fnbr[3]=4; fnbr[4]=3;
    //frot=[0,0,0,0,0];
    frot[0]=0; frot[1]=0; frot[2]=0; frot[3]=0; frot[4]=0;
    break;
  case 2:
    //fnbr=[0,5,3,1,4];
    fnbr[0]=0; fnbr[1]=5; fnbr[2]=3; fnbr[3]=1; fnbr[4]=4;
    //frot=[1,3,0,0,0];
    frot[0]=1; frot[1]=3; frot[2]=0; frot[3]=0; frot[4]=0;
    break;
  case 3:
    //fnbr=[0,5,4,2,1];
    fnbr[0]=0; fnbr[1]=5; fnbr[2]=4; fnbr[3]=2; fnbr[4]=1;
    //frot=[2,2,0,0,0];
    frot[0]=2; frot[1]=2; frot[2]=0; frot[3]=0; frot[4]=0;
    break;
  case 4:
    //fnbr=[0,5,1,3,2];
    fnbr[0]=0; fnbr[1]=5; fnbr[2]=1; fnbr[3]=3; fnbr[4]=2;
    //frot=[3,1,0,0,0];
    frot[0]=3; frot[1]=1; frot[2]=0; frot[3]=0; frot[4]=0;
    break;
  case 5:
    //fnbr=[1,3,2,4,0];
    fnbr[0]=1; fnbr[1]=3; fnbr[2]=2; fnbr[3]=4; fnbr[4]=0;
    //frot=[0,2,1,3,2];
    frot[0]=0; frot[1]=2; frot[2]=1; frot[3]=3; frot[4]=2;
    break;
  default:
    // shouldn't get here
    //fnbr=[0,0,0,0,0];
    fnbr[0]=0; fnbr[1]=0; fnbr[2]=0; fnbr[3]=0; fnbr[4]=0;
    //frot=[0,0,0,0,0];
    frot[0]=0; frot[1]=0; frot[2]=0; frot[3]=0; frot[4]=0;
    break;
  }
  return(0);
}
    
// Make sure x,y,f are within range.
// i.e. if x or y > 1 or < 0 then change face
// and recompute x,y.
// There is possible ambiguity here, so we do the best we can.
int face_range(int face, double x, double y, int *newface, double *newx, double *newy) {
  double x0,y0,tmpval; // tmp values for face, x and y
  int fnbr[5],frot[5]; // values from quadcube_arrange
  int face0,frot0; // temp rotation value

  // First check if we actually need to do anything
  if ((x >= 0)&&(x < 1)&&(y >=0)&&(y < 1)) {
    // x,y are in the proper range
    *newface=face;
    *newx=x;
    *newy=y;
    return(0);
  }  
  // set temp variables
  face0=face;
  x0=x;
  y0=y;
  //
  // do x first
  //
  quadcube_arrange(face0,fnbr,frot);
  x0=fmod(x,4);
  if (x0 < 0) {
    x0=x0+4;
  }
  if (x0 < 1) {
    // do nothing, we're already there
    frot0=0;
  } else if (x0 < 2) {
    face0=fnbr[3];
    frot0=frot[3];
    x0=x0-1;
  } else if (x0 < 3) {
    face0=fnbr[4];
    frot0=frot[4];
    x0=x0-2;
  } else {
    face0=fnbr[2];
    frot0=frot[2];
    x0=x0-3;
  }
  if (frot0 == 1) {
    tmpval=x0;
    x0=y0;
    y0=1-tmpval;
  } else if (frot0 == 2) {
    x0=1-x0;
    y0=1-y0;
  } else if (frot0 == 3) {
    tmpval=x0;
    x0=1-y0;
    y0=tmpval;
  }
  //
  // now do y
  //
  quadcube_arrange(face0,fnbr,frot);
  y0=fmod(y,4);
  if (y0 < 0) {
    y0=y0+4;
  }
  if (y0 < 1) {
    // do nothing, we're already there
    frot0=0;
  } else if (y0 < 2) {
    face0=fnbr[0];
    frot0=frot[0];
    y0=y0-1;
  } else if (y0 < 3) {
    face0=fnbr[4];
    frot0=frot[4];
    y0=y0-2;
  } else {
    face0=fnbr[1];
    frot0=frot[1];
    y0=y0-3;
  }
  if (frot0 == 1) {
    tmpval=x0;
    x0=y0;
    y0=1-tmpval;
  } else if (frot0 == 2) {
    x0=1-x0;
    y0=1-y0;
  } else if (frot0 == 3) {
    tmpval=x0;
    x0=1-y0;
    y0=tmpval;
  }
  *newface=face0;
  *newx=x0;
  *newy=y0;
  return(0);
}
