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
   squid_type nside;
   double x,y,lon,lat;
   int projection;

   if (argc != 6) {
      printf("Example usage...\n");
      printf("%s projection squid lon lat tside\n",argv[0]);
      printf("projections: 0=TSC, 1=CSC, 2=QSC, 3=HSC\n");
      printf("lon,lat in decimal degrees\n");
      printf("tside is the #pix per side of tile\n");
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
   lon=atof(argv[3]);
   lat=atof(argv[4]);
   nside=(squid_type)atoll(argv[5]);


   // Make sure squid is valid
   if (squid_validate(squid) == 0) {
      fprintf(stderr,"invalid squid argument in %s\n",argv[0]);
      exit(-1);
   }

   lon=lon*PI/180.0;
   lat=lat*PI/180.0;
   tile_sph2xy(projection,squid,lon,lat,nside,&x,&y);

   printf("x=%.5f y=%.5f\n",x,y);
   return(0);

}
