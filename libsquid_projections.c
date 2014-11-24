//
// Part of libsquid that handles the projection computations
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


//-------------------------------------------------------------------------------
// General projection calls
//
// Note: These projection algorithms are based on the following papers:
//
//       Calabretta & Greisen, "Representations of celestial coordinates in FITS",
//            A&A 395, 1077-1122, (2002)
//       Gorski, et al., "HEALPix: A Framework For High Resolution Discretization
//            and Fast Analysis of Data Distributed on the Sphere", ApJ 622, 759-771,
//            (2005)
//       Calabretta & Roukema, "Mapping on the HEALPix Grid",
//            MNRAS 381, 865-872 (2007)
//       Calabretta & Lowe, "Representing the butterfly projection in FITS --
//            projection code XPH", PASA (2013)
//
//       The Calabretta papers can all be found at
//       http://http://www.atnf.csiro.au/people/mcalabre/WCS/         
//-------------------------------------------------------------------------------

// convert face x,y to lon,lat, can handle x and y outside face limits
// nside=# pix on side of face
// x and y go from 0 to 1 across face
// lon,lat are in radians
int xyf2sph(int projection, double x, double y, int face, double *lon, double *lat) {
  double lon0; // temp variable for LON_DIR flip if needed
  if (projection == TSC) {
    if (tsc_xyf2sph(x, y, face, &lon0, lat) < 0) {
      fprintf(stderr,"tsc_xyf2sph failed in xyf2sph\n");
      return(-1);
    }
  } else if (projection == CSC) {
    if (csc_xyf2sph(x, y, face, &lon0, lat) < 0) {
      fprintf(stderr,"csc_xyf2sph failed in xyf2sph\n");
      return(-1);
    }
  } else if (projection == QSC) {
    if (qsc_xyf2sph(x, y, face, &lon0, lat) < 0) {
      fprintf(stderr,"qsc_xyf2sph failed in xyf2sph\n");
      return(-1);
    }
  } else if (projection == HSC) {
    if (hsc_xyf2sph(x, y, face, &lon0, lat) < 0) {
      fprintf(stderr,"hsc_xyf2sph failed in xyf2sph\n");
      return(-1);
    }
  } else {
    fprintf(stderr,"unknown projection in xyf2sph\n");
    return(-1);
  }
  // Flip lon if LON_DIR > 0
  // This is backwards because this library was originally written
  // for astronomy where LON_DIR = -1 is the norm.
  *lon = (LON_DIR < 0) ? lon0 : fmod(2.0*(TWOPI+LON_POLE)-lon0,TWOPI);
  return(0);
}

// convert lon,lat to face,x,y 
// x,y go from 0 to 1 across face
// lon,lat are in radians
int sph2xyf(int projection, double lon, double lat, double *x, double *y, int *face) {
  double lon0; // temp variable for LON_DIR flip if needed
  // Flip lon if LON_DIR > 0
  // This is backwards because this library was originally written
  // for astronomy where LON_DIR = -1 is the norm.
  lon0 = (LON_DIR < 0) ? lon : fmod(2.0*(TWOPI+LON_POLE)-lon,TWOPI);
  if (projection == TSC) {
    if (tsc_sph2xyf(lon0, lat, x, y, face) < 0) {
      fprintf(stderr,"tsc_sph2xyf failed in sph2xyf\n");
      return(-1);
    }
  } else if (projection == CSC) {
    if (csc_sph2xyf(lon0, lat, x, y, face) < 0) {
      fprintf(stderr,"csc_sph2xyf failed in sph2xyf\n");
      return(-1);
    }
  } else if (projection == QSC) {
    if (qsc_sph2xyf(lon0, lat, x, y, face) < 0) {
      fprintf(stderr,"qsc_sph2xyf failed in sph2xyf\n");
      return(-1);
    }
  } else if (projection == HSC) {
    if (hsc_sph2xyf(lon0, lat, x, y, face) < 0) {
      fprintf(stderr,"hsc_sph2xyf failed in sph2xyf\n");
      return(-1);
    }
  } else {
    fprintf(stderr,"unknown projection in sph2xyf\n");
    return(-1);
  }
  return(0);
}

// 
// Get number of pix on a tile side given a cdelt in deg/pix.
// k=resolution parameter
//
int squid_tside(int projection, int k, double cdelt, double *tside) { // cdelt is deg/pix
  if (projection == HSC) {
    *tside = sqrt(3/PI)*120.0/(cdelt*pow(2,k));
  } else {
    *tside = 90.0/(cdelt*pow(2,k)); // true quadcube projections
  }
  return(0);
}

//-------------------------------------------------------------------
// TSC (Tangental Spherical Cube) projection computations.
// This projection is gnomonic, but not equal area.
//-------------------------------------------------------------------

// xy2sph conversion for tangental projection
// x,y go from 0 to 1 across face.
// lon,lat are in radians
int tsc_xyf2sph(double x, double y, int face, double *lon, double *lat) {
  double X,Y; // like x,y but goes from -1 to 1 across face
  double x0, y0; // x,y but put in proper face
  int face0; // proper face
  double L,M,N; // direction cosines
  double lon0; // lon before pole_lon rotation

  // First make sure x and y are within proper limits
  if (face_range(face,x,y,&face0,&x0,&y0) < 0) {
    fprintf(stderr,"face_range failed in tsc_xyf2sph\n");
    return(-1);
  }
  // get X,Y ( -1 to 1 across face)
  X=1.0-2.0*x0;
  Y=2.0*y0-1.0;
  // get direction cosines.
  L=1/sqrt(1+X*X+Y*Y);
  M=X*L;
  N=Y*L;
  // Convert direction cosines and face to lon,lat
  LMNf2sph(L,M,N,face0,&lon0,lat);
  *lon=fmod(TWOPI+lon0+LON_POLE,TWOPI);
  return(0);
}

// sph2xy conversion for tangental projection
// x,y go from 0 to 1 across face
// lon,lat are in radians
int tsc_sph2xyf(double lon, double lat, double *x, double *y, int *face) {
  double L,M,N; // direction cosines
  double lon0; // lon before pole lon rotation
  double x0,y0; // x,y output before face range check
  int face0; // face output before range check

  lon0=fmod(TWOPI+lon-LON_POLE,TWOPI);
  // get direction cosines
  sph2LMN(lon0, lat, &face0, &L, &M, &N);
  // get x,y
  x0=(1.0-M/L)/2.0;
  y0=(N/L+1.0)/2.0;
  // Make sure x and y are within proper limits
  if (face_range(face0,x0,y0,face,x,y) < 0) {
    fprintf(stderr,"face_range failed in tsc_sph2xyf\n");
    return(-1);
  }
  return(0); // no error check needed, yet
}


//-------------------------------------------------------------------
// CSC (COBE Spherical Cube) projection computations.
// This projection is only approximately equal area,
// rms deviation ~1% over full face, 0.6% in inner 64% of face,
// but up to 14% at edge.
// The forward and backward conversions are not exact inverses,
// errors up to 24 arcsec at face edges.  
// This takes the longest to compute of the four projections.
//-------------------------------------------------------------------

// xy2sph conversion for cobe projection
// x,y go from 0 to 1 across face.
// lon,lat are in radians
int csc_xyf2sph(double x, double y, int face, double *lon, double *lat) {
  double X,Y; // like x,y but goes from -45 to 45 across face
  double x0, y0; // x,y but put in proper face
  int face0; // proper face
  double L,M,N; // direction cosines
  double lon0; // lon before pole_lon rotation
  double XW, YW; // warped versions of X,Y

  // First make sure x and y are within proper limits
  if (face_range(face,x,y,&face0,&x0,&y0) < 0) {
    fprintf(stderr,"face_range failed in tsc_xyf2sph\n");
    return(-1);
  }
  // get X,Y ( -1 to 1 across face)
  X=45.0*(1.0-2.0*x0);
  Y=45.0*(2.0*y0-1.0);
  // warp X,Y to get XW, YW
  XW=csc_rwarp(X,Y);
  YW=csc_rwarp(Y,X);
  // get direction cosines.
  L=1/sqrt(1+XW*XW+YW*YW);
  M=XW*L;
  N=YW*L;
  // Convert direction cosines and face to lon,lat
  LMNf2sph(L,M,N,face0,&lon0,lat);
  *lon=fmod(TWOPI+lon0+LON_POLE,TWOPI);
  return(0);
}

// sph2xy conversion for cobe projection
// x,y go from 0 to 1 across face
// lon,lat are in radians
int csc_sph2xyf(double lon, double lat, double *x, double *y, int *face) {
  double L,M,N; // direction cosines
  double X,Y; // like x,y but goes from -1 to 1 across face
  double XW,YW; // warped versions of X,Y
  double lon0; // lon before pole lon rotation
  double x0,y0; // output values of x,y before face range check
  int face0; // output value of face before face range check

  lon0=fmod(TWOPI+lon-LON_POLE,TWOPI);
  // get direction cosines
  sph2LMN(lon0, lat, &face0, &L, &M, &N);
  X=M/L;
  Y=N/L;
  XW=csc_fwarp(X,Y);
  YW=csc_fwarp(Y,X);
  // get x,y
  x0=(1.0-XW)/2.0;
  y0=(YW+1.0)/2.0;
  // Make sure x and y are within proper limits
  if (face_range(face0,x0,y0,face,x,y) < 0) {
    fprintf(stderr,"face_range failed in csc_sph2xyf");
    return(-1);
  }
  return(0); // no error check needed, yet
}

// CSC forward warping function
double csc_fwarp(double x, double y) {
  double g, M, G, O; // warping parameters
  double C[3][3]; // more warping parameters
  double D[2]; // yet more warping parameters
  double csum, dsum; // intermediate calculations
  int i,j; // loop counters
  double fout; // warping function output value
  
  // initialize parameters to those used by COBE project
  g=1.37484847732;
  M=0.004869491981;
  G=-0.13161671474;
  O=-0.159596235474;
  C[0][0]=0.141189631152;
  C[0][1]=-0.2815285335557;
  C[0][2]=0.106959469314;
  C[1][0]=0.0809701286525;
  C[1][1]=0.15384112876;
  C[1][2]=0.0;
  C[2][0]=-0.178251207466;
  C[2][1]=0.0;
  C[2][2]=0.0;
  D[0]=0.0759196200467;
  D[1]=-0.0217762490699;

  // compute warping factor and return value
  csum=0;
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      csum=csum+(C[i][j]*pow(x,2*i)*pow(y,2*j));
    }
  }
  csum=csum*(1-pow(y,2));
  dsum=0;
  for (i=0; i<2; i++) {
    dsum=dsum+(D[i]*pow(x,2*i));
  }
  dsum=dsum*(1-pow(x,2));
  fout = x*g + pow(x,3)*(1-g) +
         x*pow(y,2)*(1-pow(x,2))*(G+(M-G)*pow(x,2)+csum) +
         pow(x,3)*(1-pow(x,2))*(O-dsum);
  return(fout);
}

// CSC reverse warping function
double csc_rwarp(double x, double y) {
  double P[7][7]; // warping parameter matrix
  int i,j; // loop counters
  double psum; // intermediate calculation
  double xx,yy; // x/45 and y/45
  double rout; // warping function output value

  // initialize warping parameters to those used by COBE project
  for (i=0; i<7; i++) {
    for (j=0; j<7; j++) {
      P[i][j]=0.0;
    }
  }
  P[0][0]=-0.27292696;
  P[0][1]=-0.02819452;
  P[0][2]=0.27058160;
  P[0][3]=-0.60441560;
  P[0][4]=0.93412077;
  P[0][5]=-0.63915306;
  P[0][6]=0.14381585;
  P[1][0]=-0.07629969;
  P[1][1]=-0.01471565;
  P[1][2]=-0.56800938;
  P[1][3]=1.50880086;
  P[1][4]=-1.41601920;
  P[1][5]=0.52032238;
  P[2][0]=-0.22797056;
  P[2][1]=0.48051509;
  P[2][2]=0.30803317;
  P[2][3]=-0.93678576;
  P[2][4]=0.33887446;
  P[3][0]=0.54852384;
  P[3][1]=-1.74114454;
  P[3][2]=0.98938102;
  P[3][3]=0.08693841;
  P[4][0]=-0.62930065;
  P[4][1]=1.71547508;
  P[4][2]=-0.83180469;
  P[5][0]=0.25795794;
  P[5][1]=-0.53022337;
  P[6][0]=0.02584375;

  // Not computer function output
  xx=x/45.0;
  yy=y/45.0;
  psum=0.0;
  for (j=0; j<7; j++) {
    for (i=0; i<(7-j); i++) {
      psum=psum+(P[i][j]*pow(xx,2*i)*pow(yy,2*j));
    }
  }
  rout=xx+(xx*(1-pow(xx,2))*psum);
  return(rout);
}


//-------------------------------------------------------------------
// QSC (Quadrilateral Spherical Cube) projection computations.
// This projection is equal area.
//-------------------------------------------------------------------

// xy2sph conversion for quadrilateral projection
// x,y go from 0 to 1 across face.
// lon,lat are in radians
int qsc_xyf2sph(double x, double y, int face, double *lon, double *lat) {
  double X,Y; // like x,y but goes from -1 to 1 across face
  double x0, y0; // x,y but put in proper face
  int face0; // proper face
  double u,v,w; // intermediate values for dir cosine calculation
  double L,M,N; // direction cosines
  double lon0; // lon before pole_lon rotation

  // First make sure x and y are within proper limits
  if (face_range(face,x,y,&face0,&x0,&y0) < 0) {
    fprintf(stderr,"face_range failed in qsc_xyf2sph\n");
    return(-1);
  }
  // get X,Y ( -1 to 1 across face)
  X=1.0-2.0*x0;
  Y=2.0*y0-1.0;
  // determine u,v,w
  if (fabs(X) > fabs(Y)) {
    u=X; v=Y;
  } else {
    u=Y; v=X;
  }
  w = (u == 0) ? 0.0 : sin((PI/12.0)*v/u)/(cos((PI/12.0)*v/u)-1/sqrt(2));
  // get direction cosines.
  L=1.0-(pow(u,2)*(1.0-1.0/sqrt(2.0+pow(w,2))));
  if (u == X) {
    M=sqrt((1.0-pow(L,2))/(1.0+pow(w,2)));
    if (X < 0) M=-1*M;
    N=M*w;
  } else {
    N=sqrt((1.0-pow(L,2))/(1.0+pow(w,2)));
    if (Y < 0) N=-1*N;
    M=N*w;
  }
  // Convert direction cosines and face to lon,lat
  LMNf2sph(L,M,N,face0,&lon0,lat);
  *lon=fmod(TWOPI+lon0+LON_POLE,TWOPI);
  return(0);
}

// sph2xy conversion for quadrilateral projection
// x,y go from 0 to 1 across face
// lon,lat are in radians
int qsc_sph2xyf(double lon, double lat, double *x, double *y, int *face) {
  double L,M,N; // direction cosines
  double u,v,w,s; // intermediate values for calculation
  double lon0; // lon before pole lon rotation
  double x0,y0; // x,y output before face range check
  int face0; // face output before range check

  lon0=fmod(TWOPI+lon-LON_POLE,TWOPI);
  // get direction cosines
  sph2LMN(lon0, lat, &face0, &L, &M, &N);
  // get u,v,w,s
  s = ((M > fabs(N))||(N >= fabs(M))) ? 1.0 : -1.0;
  if (fabs(M) > fabs(N)) {
    if (M == 0) {
	u = s*sqrt(1.0-L);
	v = (N >= 0) ? 3.0*u : -3.0*u;
    } else {
      w = N/M;
      u = s*sqrt((1.0-L)/(1.0-1.0/sqrt(2.0+pow(w,2))));
      v = u*(atan(w)-asin(w/sqrt(2.0*(1.0+pow(w,2)))))/(DD2R*15.0);
    }
  } else {
    if (N == 0) {
	u = s*sqrt(1.0-L);
	v = (M >= 0) ? 3.0*u : -3.0*u;
    } else {
      w = M/N;
      u = s*sqrt((1.0-L)/(1.0-1.0/sqrt(2.0+pow(w,2))));
      v = u*(atan(w)-asin(w/sqrt(2.0*(1.0+pow(w,2)))))/(DD2R*15.0);
    }
  } 
  // get x,y
  if (fabs(M) > fabs(N)) {
    x0=1.0-(1.0+u)/2.0;
    y0=(1.0+v)/2.0;
  } else {
    x0=1.0-(1.0+v)/2.0;
    y0=(1.0+u)/2.0;
  }
  // Make sure x and y are within proper limits
  if (face_range(face0,x0,y0,face,x,y) < 0) {
    fprintf(stderr,"face_range failed in qsc_sph2xyf\n");
    return(-1);
  }
  return(0); // no error check needed, yet
}


//-------------------------------------------------------------------
// HSC (HEALPix Spherical Cube) projection computations.
// This projection is equal area, but actually a hybrid projection.
// The equatorial faces are cylindrical equal area.
// The polar faces are Collignon projections.
// In general, this is the fastest to calculate of the 4 projections,
// often by a lot.
//-------------------------------------------------------------------

// xy2sph conversion, can handle x and y outside face limits
// x,y go from 0 to 1 across face
// lon,lat are in radians
int hsc_xyf2sph(double x, double y, int face, double *lon, double *lat) {
  double x1,y1,x2,y2;
  double alpha, alpha1, r1;
  double x0,y0; // x,y put on proper face
  int face0; // proper face

  // First make sure x and y are within proper limits
  if (face_range(face,x,y,&face0,&x0,&y0) < 0) {
    fprintf(stderr,"face_range failed in xyf2sph\n");
    return(-1);
  }
  face=face0;
  x=x0;
  y=y0;

  *lon=-1.0; // error value
  *lat=-1.0; // error value
  if ((face == 0)||(face == 5)) {
    // Polar region
    x2=1.0-2.0*x;
    y2=1.0-2.0*y;
    alpha=atan2(x2,y2); // yes, x and y are reversed
    alpha1=fmod(alpha+QUARTPI+(2*PI),HALFPI)-QUARTPI;
    r1=sqrt(x2*x2+y2*y2);
    if (face == 0) {
      // North polar region
      if (r1 == 0) {
	*lon=0.0;
	*lat=HALFPI;
      } else {
	x1=r1*sin(alpha1);
	y1=r1*cos(alpha1); // also = sigma
	*lon=(alpha-alpha1)+(x1/y1)*QUARTPI+LON_POLE; // center of face1 is defined by LON_POLE
	*lat=asin(1.0-((y1*y1)/3.0));
      }
    } else {
      // South polar region
      if (r1 == 0) {
	*lon=0.0;
	*lat=-1.0*HALFPI;
      } else {
	x1=r1*sin(alpha1);
	y1=-1.0*r1*cos(alpha1); // also = sigma
	*lon=(PI-alpha+alpha1)+(x1/y1)*QUARTPI+LON_POLE; // center of face1 is defined by LON_POLE
	*lat=-1.0*asin(1.0-((y1*y1)/3.0));
      }
    }
  } else {
    // Equatorial region
    *lon=HALFPI*(1.0*face-x)-QUARTPI+LON_POLE; // center of face1 is defined by LON_POLE
    y1=-1.0*QUARTPI+(HALFPI*y);
    *lat=asin(y1/(HALFPI*0.75));
  }
  *lon=fmod(*lon+2*PI,2*PI);

  return(0);
}

//    
// sph2xy conversion for healpix projection
// x,y go from 0 to 1 across face
// lon,lat are in radians
// 
int hsc_sph2xyf(double lon, double lat, double *x, double *y, int *face) {
  double x0, y0, face0, x1, y1, lon1, lat1, r1, alpha, alpha1, face1;
  
  if (fabs(lat) > THETAX) {
    // Polar region
    y1=sqrt(3.0*(1.0-fabs(sin(lat)))); // also = sigma
    lon1=fmod(lon+QUARTPI+(2*PI)+LON_POLE,HALFPI)-QUARTPI; // center of face1 is defined by LON_POLE
    x1=(lon1*y1)/QUARTPI;
    r1=sqrt(x1*x1+y1*y1);
    alpha1=atan2(x1,y1); // yes, x and y are reversed
    alpha=alpha1+(lon-lon1)-LON_POLE; // center of face1 is defined by LON_POLE
    if (lat > 0) {
      // North polar region
      face0=0;
      x0=0.5-(r1*sin(alpha)/2.0);
      y0=0.5-(r1*cos(alpha)/2.0);
    } else {
      // South polar region
      face0=5;
      x0=0.5-(r1*sin(alpha)/2.0);
      y0=0.5+(r1*cos(alpha)/2.0);
    }
  } else {
    // Equatorial region
    lon1=fmod(lon+QUARTPI+(2*PI)-LON_POLE,2*PI); // center of face1 is defined by LON_POLE
    lat1=HALFPI*sin(lat)*3.0/4.0+QUARTPI;
    face1=floor(lon1/HALFPI)+1;
    x0=1.0*face1-lon1/HALFPI;
    y0=lat1/HALFPI;
    face0=(int)face1;
  }
  // Make sure x and y are within proper limits
  if (face_range(face0,x0,y0,face,x,y) < 0) {
    fprintf(stderr,"face_range failed in hsc_sph2xyf\n");
    return(-1);
  }
  return(0);
}

//-------------------------------------------------------------------
// General projection utility functions
//-------------------------------------------------------------------

// Get the face number for a given lon,lat for TSC, CSC, and QSC projections.
// Note that this doesn't work for the HSC projection!
// lon,lat in londians
void quadcube_getface(double lon, double lat, int *face) {
  int i, tmpface=0; // loop counter, temp face value
  double sdist, mindist=1000; // sph dist to face centers
  double lonf; // lon of face
  // do north pole
  sphdist(lon,lat,LON_POLE,HALFPI,&sdist);
  if (sdist < mindist) {
    mindist=sdist;
    tmpface=0;
  }
  // do south pole
  sphdist(lon,lat,LON_POLE,-1.0*HALFPI,&sdist);
  if (sdist < mindist) {
    mindist=sdist;
    tmpface=5;
  }
  // do equator
  for (i=1; i<5; i++) {
    lonf=fmod(LON_POLE+((i-1)*HALFPI),TWOPI);
    sphdist(lon,lat,lonf,0.0,&sdist);
    if (sdist < mindist) {
      mindist=sdist;
      tmpface=i;
    }
  }
  *face=tmpface;
}

// Get face number for a given lon,lat for HSC projection.
// lon,lat are in radians
void hsc_getface(double lon, double lat, int *face) {
  double londiff; // difference between lon and LON_POLE + 45deg
  lon=fmod(TWOPI+lon,TWOPI); // make sure lon in proper range
  if (lat > THETAX) {
    *face=0;
  } else if (lat < (-1*THETAX)) {
    *face=5;
  } else {
    londiff=fmod(TWOPI+lon-LON_POLE+QUARTPI, TWOPI);
    if (londiff <= HALFPI) {
      *face=1;
    } else if (londiff <= PI) {
      *face=2;
    } else if (londiff <= THREEHALFPI) {
      *face=3;
    } else {
      *face=4;
    }
  }
}

//
// Get direction cosines for a particular face.
// face is calculated based on lon,lat and returned.
// lon,lat are in radians
// L,M,N are the direction cosines for the face
// corresponding to lon,lat
//
int sph2LMN(double lon, double lat, int *face, double *L, double *M, double *N) {
  double l,m,n; // direction cosines

  // get face
  quadcube_getface(lon,lat,face);

  // get direction cosines
  l=cos(lat)*cos(lon);
  m=cos(lat)*sin(lon);
  n=sin(lat);

  // get direction cosines for this face
  switch (*face) {
  case 0:
    *M=m;
    *N=-1*l;
    *L=n;
    break;
  case 1:
    *M=m;
    *N=n;
    *L=l;
    break;
  case 2:
    *M=-1*l;
    *N=n;
    *L=m;
    break;
  case 3:
    *M=-1*m;
    *N=n;
    *L=-1*l;
    break;
  case 4:
    *M=l;
    *N=n;
    *L=-1*m;
    break;
  case 5:
    *M=m;
    *N=l;
    *L=-1*n;
    break;
  default: // use face 1 params
    *M=m;
    *N=n;
    *L=l;
  }
  return(0);
}

//
// Convert direction cosines with face to lon,lat in radians
//
int LMNf2sph(double L, double M, double N, int face, double *lon, double *lat) {
  double l,m,n; // general direction cosines

  switch (face) {
  case 0:
    l=-1*N;
    m=M;
    n=L;
    break;
  case 1:
    l=L;
    m=M;
    n=N;
    break;
  case 2:
    l=-1*M;
    m=L;
    n=N;
    break;
  case 3:
    l=-1*L;
    m=-1*M;
    n=N;
    break;
  case 4:
    l=M;
    m=-1*L;
    n=N;
    break;
  case 5:
    l=N;
    m=M;
    n=-1*L;
    break;
  default:  // use face 1 params
    l=L;
    m=M;
    n=N;
  }
  *lat=asin(n);
  *lon=atan2(m,l);
  return(0); // no error check, yet
}


