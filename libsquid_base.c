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

#include <libsquid.h>

//
// Is this a vaild squid? Return 1 if so, 0 if not
//
int squid_validate(squid_type squid) {
  // Must be at least 4 bits long and an even number of bits
  long nbits, nshift;
  int face;
  nbits=(long)(floor(log(squid)/log(2)))+1;
  if (nbits % 2) { // number of bits not even
    return(0);
  } else if (nbits < 4) { // squid too small
    return(0);
  }
  nshift=nbits-4;
  if (nshift > 0) {
    face=(squid >> nshift)-8;
  } else {
    face=squid-8;
  }
  if (face > 5) { // invalid face
    return(0);
  }
  return(1);
}

//
// Get squid resolution parameter.
// Returns a -1 if not a valid squid.
//
int squid_getres(squid_type squid) {
  int k; // resolution parameter
  double k1; // temporary double 
  if (squid_validate(squid) == 0) { // squid not valid!
    return(-1);
  }
  k1=floor((log(squid)/log(4))-1);
  k=(int)k1;
  return(k);
}

//
// Get the face number corresponding to a squid
// Returns -1 if not a valid squid
//
int squid_getface(squid_type squid) {
  int face, k,retval;
  squid_type x,y;
  if (squid_validate(squid) == 0) { // squid not valid!
    fprintf(stderr,"invalid squid in squid_getface\n");
    return(-1);
  }
  retval=squid2xyfk(squid,&x,&y,&face,&k);
  if (retval >= 0) {
    return(face);
  } else {
    fprintf(stderr,"squid2xyfk failed in squid_getface\n");
    return(-1);
  }
}

//
// Convert squid to x,y,f,k where:
// f=face, k=resolution parameter
// x,y = quad-tree coords within face
//
int squid2xyfk(squid_type squid, squid_type *x, squid_type *y, int *face, int *k) {
  int k1;
  squid_type mask, qtid;
  if (squid_validate(squid) == 0) {
    fprintf(stderr,"invalid squid in squid2xyfk\n");
    return(-1);
  }
  k1=squid_getres(squid);
  mask=(squid_type)(pow(2,2*k1)-1);
  qtid = squid & mask;
  *face = (int)((squid | mask)>>(2*k1))-8;
  *k = k1;
  if (qtree_id2xy(qtid, k1, x, y) == -1) {
    fprintf(stderr,"qtree_id2xy failed in squid2xyfk\n");
    return(-1);
  }
  return(0);
}

//
// Convert x,y,f,k to a squid where:
// f=face, k=resolution parameter
// x,y = quad-tree coords within face
//
int xyfk2squid(squid_type x, squid_type y, int face, int k, squid_type *squid) {
  squid_type qtid, face1;
  qtree_xy2id(x,y,k,&qtid);
  if (qtid == -1) {
    fprintf(stderr,"qtree_xy2id failed in xyfk2squid\n");
    return(-1);
  }
  face1 = (squid_type)face + 8;
  *squid=(face1 << (2*k))+qtid;
  return(0);
}
  
//  
// Get lon,lat of squid center.
// lon,lat in radians
//
int squid2sph(int projection, squid_type squid, double *lon, double *lat) {
  int k, face;
  squid_type x,y;
  double x1,y1;

  if (squid2xyfk(squid, &x, &y, &face, &k) == -1) {
    fprintf(stderr,"squid2xyfk failed in squid2sph\n");
    return(-1);
  }
  x1=((double)x+0.5)/pow(2,k); // convert to face x
  y1=((double)y+0.5)/pow(2,k); // convert to face y
  if (xyf2sph(projection,x1,y1,face,lon,lat) == -1) {
    fprintf(stderr,"xyf2sph failed in squid2sph\n");
    return(-1);
  }
  return(0);
}

//
// Get squid corresponding to lon,lat,k where:
// lon,lat in radians, k=resolution parameter
//
int sph2squid(int projection, double lon, double lat, int k, squid_type *squid) {
  int face;
  squid_type xl,yl; // quad-tree coords
  double x,y;

  // convert lon,lat to face x,y
  if (sph2xyf(projection,lon,lat,&x,&y,&face) == -1) {
    fprintf(stderr,"sph2xyf failed in sph2xyf\n");
    return(-1);
  }
  if (x == 1.0) {
    // make sure we stay on proper face
    x=x-1e-10;
  }
  if (y == 1.0) {
    // make sure we stay on proper face
    y=y-1e-10;
  }
  xl=(squid_type)floor(x*pow(2,k)); // quad-tree coord for this resolution
  yl=(squid_type)floor(y*pow(2,k)); // quad-tree coord for this resolution
  if (xyfk2squid(xl,yl,face,k,squid) == -1) {
    fprintf(stderr,"xyfk2squid failed in sph2squid\n");
    return(-1);
  }
  return(0);
}

//
// For a given squid tile and x,y,tside get lon,lat
//
int tile_xy2sph(int projection, squid_type squid, double x, double y, squid_type tside, double *lon, double *lat) {
  squid_type x0,y0;
  double x1,y1;
  int k,face;
  
  if (squid2xyfk(squid,&x0,&y0,&face,&k) == -1) {
    fprintf(stderr,"squid2xyfk failed in tile_xy2sph\n");
    return(-1);
  }
  x1=(x/tside+x0)/pow(2,k); // convert to face x
  y1=(y/tside+y0)/pow(2,k); // convert to face y
  if (x1 == 1.0) {
    // make sure we stay on proper face
    x1=x1-1e-10;
  }
  if (y1 == 1.0) {
    // make sure we stay on proper face
    y1=y1-1e-10;
  }
  if (xyf2sph(projection,x1,y1,face,lon,lat) == -1) {
    fprintf(stderr,"xyf2sph failed in tile_xy2sph\n");
    return(-1);
  }
  return(0);
}

//
// For a given squid tile and lon,lat,tside get x,y
//
int tile_sph2xy(int projection, squid_type squid, double lon, double lat, squid_type tside, double *x, double *y) {
  squid_type x0,y0;
  double x1,y1;
  squid_type nside1;
  int k,face0,face1;

  if (squid2xyfk(squid,&x0,&y0,&face0,&k) == -1) {
    fprintf(stderr,"squid2xyfk failed in tile_sph2xy\n");
    return(-1);
  }
  nside1=(squid_type)pow(2,k)*tside;

  if (sph2xyf(projection,lon,lat,&x1,&y1,&face1) == -1) {
    fprintf(stderr,"sph2xyf failed in tile_sph2xy\n");
    return(-1);
  }
  *x=nside1*x1-tside*x0;
  *y=nside1*y1-tside*y0;

  return(0);
}

//
// get 4 element lon and lat arrays with corner coordinates for squid tile
//
int squid_corners(int projection, squid_type squid, double *lonrr, double *darr) {
  double nside;

  nside=1000.0;
  if (tile_xy2sph(projection,squid,0.0,0.0,(squid_type)nside,&lonrr[0],&darr[0]) == -1) {
    fprintf(stderr,"tile_xy2sph failed in squid_corners\n");
    return(-1);
  }
  if (tile_xy2sph(projection,squid,nside,0.0,(squid_type)nside,&lonrr[1],&darr[1]) == -1) {
    fprintf(stderr,"tile_xy2sph failed in squid_corners\n");
    return(-1);
  }
  if (tile_xy2sph(projection,squid,0.0,nside,(squid_type)nside,&lonrr[2],&darr[2]) == -1) {
    fprintf(stderr,"tile_xy2sph failed in squid_corners\n");
    return(-1);
  }
  if (tile_xy2sph(projection,squid,nside,nside,(squid_type)nside,&lonrr[3],&darr[3]) == -1) {
    fprintf(stderr,"tile_xy2sph failed in squid_corners\n");
    return(-1);
  }

  return(0);
}

// Get nearest point on tile border to sky coords
// Squid is id of tile
// lons, lats are coords of search point in rad
// lonn, latn are coords of nearest point on tile border in rad
int tile_nearest(int projection, squid_type squid, double lons, double lats, double *lonn, double *latn) {
  squid_type nside;
  double tinyval, nsideh;
  double lonc[4], latc[4]; // squid corner coords, corner distances
  double cdmin; // minimum corner distance
  long i, j, i0, j0, i1; // loop counters, index of closest corner and midpoint
  double sdist; // spherical distance in rad
  double midx[4],midy[4]; // x,y of midpionts
  double mdmin; // minimum dist to midpoints
  double lonm[4],latm[4],mdarr[4]; // lon,lat of midpoints, array of midpoint distances
  double lonm0, latm0; // temp midpoint lon,lat
  double rc0,dc0,rc1,dc1,rm0,dm0; // cords of closest corners and midpoint
  double xt,yt; // xy tile coords of midpoint
  double lonn0, latn0; // temp values of lonn,latn
  double xt0,yt0; // tile xy of nearest point
  int maxiter; // max number of iterations
  double width, rcutoff, lon0chk; // for checks within iteration loop
  double xt1, yt1, xt1_tmp, yt1_tmp; // temp midpoint values within iteration loop
  double lonni, latni;  // values of lonn,latn before iteration starts
  double lonn1, latn1; // temp values for lonn,latn within iteration loop
  double sdist1; // sph dist between iteration solutions
  double initial_err; // error (sdist) from initial solution to iterated solution
  
  tinyval=1e-6;
  nside=1000;
  nsideh=nside/2.0;

  // get corners
  if (squid_corners(projection, squid, lonc, latc) < 0) {
    fprintf(stderr,"squid_corners failed in tile_nearest\n");
    return(-1);
  }
  // find closest corner
  cdmin=1000.0;
  i0=0;
  for (i=0; i<4; i++) {
    sphdist(lonc[i],latc[i],lons,lats,&sdist);
    if (sdist < cdmin) {
      cdmin=sdist;
      i0=i;
    }
  }
  // find closest midpoint
  midx[0]=nsideh; midx[1]=0.0; midx[2]=nside; midx[3]=nsideh;
  midy[0]=0.0; midy[1]=nsideh; midy[2]=nsideh; midy[3]=nside;
  mdmin=1000.0;
  j0=0;
  for (j=0; j<4; j++) {
    if (tile_xy2sph(projection, squid, midx[j], midy[j], nside, &lonm0, &latm0) < 0) {
      fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
      return(-1);
    }
    lonm[j]=lonm0;
    latm[j]=latm0;
    sphdist(lonm0,latm0,lons,lats,&sdist);
    mdarr[j]=sdist;
    if (sdist < mdmin) {
      mdmin=sdist;
      j0=j;
    }
  }

  // Get best midpoint and corners
  if (cdmin <= mdmin) { // corner is closer, get best midpoint and opposite corner
    if (i0 == 0) {
      if (mdarr[0] < mdarr[1]) {
        j0=0;
        i1=1;
      } else {
        j0=1;
        i1=2;
      }
    } else if (i0 == 1) {
      if (mdarr[0] < mdarr[2]) {
        j0=0;
        i1=0;
      } else {
        j0=2;
        i1=3;
      }
    } else if (i0 == 2) {
      if (mdarr[1] < mdarr[3]) {
        j0=1;
        i1=0;
      } else {
	j0=3;
	i1=3;
      }
    } else {
      if (mdarr[2] < mdarr[3]) {
	j0=2;
	i1=1;
      } else {
	j0=3;
	i1=2;
      }
    }
  } else { // midpoint is closer, get two adjacent corners
    if (j0 == 0) {
      i0=0;
      i1=1;
    } else if (j0 == 1) {
      i0=0;
      i1=2;
    } else if (j0 == 2) {
      i0=1;
      i1=3;
    } else {
      i0=2;
      i1=3;
    }
  }
  rc0=lonc[i0];
  dc0=latc[i0];
  rc1=lonc[i1];
  dc1=latc[i1];
  rm0=lonm[j0];
  dm0=latm[j0];
  xt=midx[j0];
  yt=midy[j0];
  
  // get arc intercept
  arc_intercept(rc0, dc0, rm0, dm0, rc1, dc1, lons, lats, &lonn0, &latn0);
  if (tile_sph2xy(projection, squid, lonn0, latn0, nside, &xt0, &yt0) < 0) {
    fprintf(stderr,"tile_sph2xy failed in tile_nearest\n");
    return(-1);
  }
  //printf("Initial solution lonn=%.5f latn=%.5f\n",lonn0,latn0);

  // Now do iteration.  If arc points lie on a Cartesian plane, then
  // the above solution is exact and the iteration should exit after
  // a single loop.  
  lonni=lonn0; // preserve initial solution
  latni=latn0; // preserve initial solution
  maxiter=10;
  width=1.0; // width of arc intercept recalculation
  rcutoff=1e-6; // loop cutoff value in radians (1e-6 ~ .02")
  lon0chk=6.0/nside; // check if we are really close to diagonal, problem for HSC
  if ((fabs(lonn0) < lon0chk)||(fabs(lonn0-PI) < lon0chk)||(fabs(lonn0-TWOPI) < lon0chk)) {
    // close to diagonal, change params to finger grain
    maxiter=50;
    width=0.1;
    rcutoff=5e-6;
  }
  for (i=0; i < maxiter; i++) {
    if (xt0 < tinyval) {
      xt0=xt0+1;
    } else if (xt0 > (nside-tinyval)) {
      xt0=xt0-1;
    }
    if (yt0 < tinyval) {
      yt0=yt0+1;
    } else if (yt0 > (nside-tinyval)) {
      yt0=yt0-1;
    }
    if ((j0 == 0)||(j0 == 3)) { // either top or bottom side
      if (tile_xy2sph(projection, squid, xt0-width, yt, nside, &rc0, &dc0) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
      if (tile_xy2sph(projection, squid, xt0+width, yt, nside, &rc1, &dc1) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
    } else { // either left or right sides
      if (tile_xy2sph(projection, squid, xt, yt0-width, nside, &rc0, &dc0) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
      if (tile_xy2sph(projection, squid, xt, yt0+width, nside, &rc1, &dc1) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
    }
    if ((lonn0 == rc0)||(lonn0 == rc1)) { // make sure we get 3 distinct points
      lonn0=(rc0+rc1)/2.0;
    }
    // Re-calculate intercept with new mid and end points
    arc_intercept(rc0, dc0, lonn0, latn0, rc1, dc1, lons, lats, &lonn1, &latn1);
    if (tile_sph2xy(projection, squid, lonn1, latn1, nside, &xt1_tmp, &yt1_tmp) < 0) {
      fprintf(stderr,"tile_sph2xy failed in tile_nearest\n");
      return(-1);
    }
    // Adjust answer to make sure output is right on the line
    if ((j0 == 0)||(j0 == 3)) { // either bottom or top side
      if (tile_xy2sph(projection, squid, xt1_tmp, yt, nside, &lonn1, &latn1) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
    } else { // either left or right side
      if (tile_xy2sph(projection, squid, xt, yt1_tmp, nside, &lonn1, &latn1) < 0) {
	fprintf(stderr,"tile_xy2sph failed in tile_nearest\n");
	return(-1);
      }
    }
    // Run again with adjusted answer
    if (tile_sph2xy(projection, squid, lonn1, latn1, nside, &xt1, &yt1) < 0) {
      fprintf(stderr,"tile_sph2xy failed in tile_nearest\n");
      return(-1);
    }
    sphdist(lonn1,latn1,lonn0,latn0,&sdist1); // distance from previous solution
    //printf("sdist1=%.30f arcseconds\n",3600.0*sdist1/DD2R);
    if (sdist1 < rcutoff) { // we're done!
      break;
    }
    // update loop values for next iteration
    lonn0=lonn1;
    latn0=latn1;
    xt0=xt1;
    yt0=yt1;
  }
  sphdist(lonni,latni,lonn0,latn0,&initial_err);
  //printf("initial error=%.30f arcseconds\n",3600.0*initial_err/DD2R);
  //printf("lonn=%.5f latn=%.5f\n",lonn0/DD2R,latn0/DD2R);
  *lonn=lonn0;
  *latn=latn0;
  return(0);
}


// Get list of squid tiles within a cone region.
// lon,lat = search point in radians
// srad = search radius in radians
// kmin, kmax = the min and max squid resolutions for the search.
// Return values are:
// nfull = the number of elements in full_tiles[]
// full_tiles = pointer to array of squids that are entirely within cone
// npart = the number of elements in part_tiles[]
// part_tiles = pointer to array of squids that are partially within cone
// NOTE: You will need to free the full_tiles and part_tiles arrays eventually.
int cone_search(int projection, double lon, double lat, double srad, int kmin, int kmax,
		long *nfull, squid_type **full_tiles, long *npart, squid_type **part_tiles) {
  int i,j,k; // loop counters
  int ci; // corner index
  int ktmp; // temporary value for resolution param
  squid_type squid0, squid1, squid11; // squids corresponding to lon,lat at various resolutions
  squid_type squid_tmp; // candidate squid to add to arrays
  double lonc[4], latc[4]; // squid corner arrays
  int ccount; // how many corners are within cone
  double lonn,latn; // nearest point on tile to lon,lat
  double sdist; // spherical distance in rad
  squid_type *full_tiles0, *part_tiles0, *part_tiles1;  // temporary tile arrays
  long nfull0=0, npart0=0, npart1=0; // number of elements in respective temp arrays
  
  if (kmin > kmax) { // make sure args are in right order
    ktmp=kmax;
    kmax=kmin;
    kmin=ktmp;
  }
  if (kmax > KMAX) { // stay within squid resolution limits
    kmax=KMAX;
  }

  // Initialize pointers
  full_tiles0=malloc(sizeof(squid_type));
  if (full_tiles0 == NULL) {
    fprintf(stderr,"malloc failed in cone_search, %s\n",strerror(errno));
    return(-1);
  }
  part_tiles0=malloc(sizeof(squid_type));
  if (part_tiles0 == NULL) {
    fprintf(stderr,"malloc failed in cone_search, %s\n",strerror(errno));
    return(-1);
  }
  
  //
  // First, check the 6 face tiles
  //
  if (sph2squid(projection, lon, lat, 0, &squid0) < 0) {
    fprintf(stderr,"sph2squid failed in cone_search\n");
    return(-1);
  }
  //printf("squid0=%ld\n",squid0);
  for (squid_tmp=8; squid_tmp < 14; squid_tmp++) {
    // First check if all 4 corners are in search radius.
    // If so, it should be a fully contained tile.
    if (squid_corners(projection, squid_tmp, lonc, latc) < 0) {
      fprintf(stderr,"squid_corners failed in cone_search\n");
      return(-1);
    }
    ccount=0;
    for (ci=0; ci < 4; ci++) {
      sphdist(lonc[ci],latc[ci],lon,lat,&sdist);
      if (sdist < srad) {
	ccount++;
      }
    }
    if (ccount == 4) { // all 4 points are in!
      if (kmin == 0) {
	// we are within resolution range so append to full_tiles[]
	nfull0++;
	full_tiles0=realloc(full_tiles0, nfull0*sizeof(squid_type));
	if (full_tiles0 == NULL) {
	  fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	  return(-1);
	}
	full_tiles0[nfull0-1]=squid_tmp;
	continue;
      } else {
	// we are outside resolution range so append to part_tiles[]
	npart0++;
	part_tiles0=realloc(part_tiles0, npart0*sizeof(squid_type));
	if (part_tiles0 == NULL) {
	  fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	  return(-1);
	}
	part_tiles0[npart0-1]=squid_tmp;
	continue;
      }
    }
    // Next check if search point is within tile
    if (squid_tmp == squid0) {
      // search region within tile, append to part_tiles[]
      npart0++;
      part_tiles0=realloc(part_tiles0, npart0*sizeof(squid_type));
      if (part_tiles0 == NULL) {
	fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	return(-1);
      }
      part_tiles0[npart0-1]=squid_tmp;
      continue;
    }      
    // Finally check if search region overlaps tile
    if (tile_nearest(projection, squid_tmp, lon, lat, &lonn, &latn) < 0) {
      fprintf(stderr,"tile_nearest failed in cone_search\n");
      return(-1);
    }
    sphdist(lon,lat,lonn,latn,&sdist);
    if (sdist < srad) {
      // search region overlaps tile, append to part_tiles[]
      npart0++;
      part_tiles0=realloc(part_tiles0, npart0*sizeof(squid_type));
      if (part_tiles0 == NULL) {
	fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	return(-1);
      }
      part_tiles0[npart0-1]=squid_tmp;
    }
  }
  // We definitely should have at least one part_tile in the queue

  //
  // Loop over the resolution levels
  //
  for (k=0; k<kmax; k++) {
    //printf("Checking level %d\n",k+1);
    if (sph2squid(projection, lon, lat, k+1, &squid0) < 0) {
      fprintf(stderr,"sph2squid failed in cone_search\n");
      return(-1);
    }
    npart1=0;
    part_tiles1=malloc(sizeof(squid_type)); // temp array
    if (part_tiles1 == NULL) {
      fprintf(stderr,"malloc failed in cone_search, %s\n",strerror(errno));
      return(-1);
    }
    // Now drop down a level and check subtiles of
    // each partial squid
    for (i=0; i<npart0; i++) { // loop over existing partial tile array
      squid_tmp=part_tiles0[i];
      squid1 = squid_tmp << 2; // shift over to go down resolution level
      for (j=0; j<4; j++) {
	squid11=squid1+j; // get sub squid number
	// First check if all 4 corners are in search radius.
	// If so, it should be a fully contained tile
	if (squid_corners(projection, squid11, lonc, latc) < 0) {
	  fprintf(stderr,"squid_corners failed in cone_search\n");
	  return(-1);
	}
	ccount=0;
	for (ci=0; ci < 4; ci++) {
	  sphdist(lonc[ci],latc[ci],lon,lat,&sdist);
	  if (sdist < srad) {
	    ccount++;
	  }
	}
	if (ccount == 4) { // all 4 points are in!
	  if (kmin <= k+1) {
	    // we are within resolution range so append to full_tiles[]
	    nfull0++;
	    full_tiles0=realloc(full_tiles0, nfull0*sizeof(squid_type));
	    if (full_tiles0 == NULL) {
	      fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	      return(-1);
	    }
	    full_tiles0[nfull0-1]=squid11;
	    continue;
	  } else {
	    // we are outside resolution range so append to part_tiles[]
	    npart1++;
	    part_tiles1=realloc(part_tiles1, npart1*sizeof(squid_type));
	    if (part_tiles1 == NULL) {
	      fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	      return(-1);
	    }
	    part_tiles1[npart1-1]=squid11;
	    continue;
	  }
	}
	// Next check if search point is within tile
	if (squid11 == squid0) {
	  // append to part_tiles[]
	  npart1++;
	  part_tiles1=realloc(part_tiles1, npart1*sizeof(squid_type));
	  if (part_tiles1 == NULL) {
	    fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	    return(-1);
	  }
	  part_tiles1[npart1-1]=squid11;
	  continue;
	}
	// Finally check if search region overlaps tile
	if (tile_nearest(projection, squid11, lon, lat, &lonn, &latn) < 0) {
	  fprintf(stderr,"tile_nearest failed in cone_search, %s\n",strerror(errno));
	  return(-1);
	}
	sphdist(lon,lat,lonn,latn,&sdist);
	if (sdist < srad) {
	  // search region overlaps tile, append to part_tiles[]
	  npart1++;
	  part_tiles1=realloc(part_tiles1, npart1*sizeof(squid_type));
	  if (part_tiles1 == NULL) {
	    fprintf(stderr,"realloc failed in cone_search, %s\n",strerror(errno));
	    return(-1);
	  }
	  part_tiles1[npart1-1]=squid11;
	}
      }
    }
    // Now update part_tiles[]
    if (part_tiles0 != NULL) {
      free(part_tiles0);
    }
    npart0=npart1;
    part_tiles0=part_tiles1;
  }

  // update output pointers
  *nfull=nfull0;
  *full_tiles=full_tiles0;
  *npart=npart0;
  *part_tiles=part_tiles0;

  return(0);

}
