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

#define _GNU_SOURCE 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

#include <libsquid.h>

int main(int argc, char *argv[]) {
  squid_type squid;
  squid_type xl,yl;
  double lon,lat;
  int i,k,face;
  double lonc[4],latc[4];
  int projection;

  if (argc != 3) {
    printf("Example usage...\n");
    printf("%s projection squid\n",argv[0]);
    printf("projections: 0=TSC, 1=CSC, 2=QSC, 3=HSC\n");
    exit(-1);
  }
  projection=atoi(argv[1]);
  if (projection == TSC) {
    printf("TSC Projection\n");
  } else if (projection == CSC) {
    printf("CSC Projection\n");
  } else if (projection == QSC) {
    printf("QSC Projection\n");
  } else if (projection == HSC) {
    printf("HSC Projection\n");
  } else {
    printf("Unknown projection! Using HSC.\n");
    projection=HSC;
  }
  squid=(squid_type)atoll(argv[2]);

  // Make sure squid is valid
  if (squid_validate(squid) == 0) {
    fprintf(stderr,"invalid squid argument in %s\n",argv[0]);
    exit(-1);
  }

  // Get basic squid information
  if (squid2xyfk(squid,&xl,&yl,&face,&k) == -1) {
    fprintf(stderr,"squid2xyfk failed in %s\n",argv[0]);
    exit(-1);
  }

  // Get squid center
  if (squid2sph(projection,squid,&lon,&lat) == -1) {
    fprintf(stderr,"squid2sph failed in %s\n",argv[0]);
    exit(-1);
  }

  // Get squid corners
  if (squid_corners(projection,squid,lonc,latc) == -1) {
    fprintf(stderr,"squidcorners failed in %s\n",argv[0]);
    exit(-1);
  }

  // Print out information
  printf("SQUID=%ld\n",(long)squid);
  printf("LEVEL=%d, FACE=%d, X=%ld, Y=%ld\n",k,face,(long)xl,(long)yl);
  printf("CENTER LON=%.6f LAT=%.6f\n",lon/DD2R,lat/DD2R);
  for (i=0; i<4; i++) {
    printf("CORNER %d LON=%.6f LAT=%.6f\n",i,lonc[i]/DD2R,latc[i]/DD2R);
  }
    
  return(0);

}
