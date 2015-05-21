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
   double lon,lat,rad;
   long nfull, npart;
   squid_type *full_tiles, *part_tiles;
   long i;
   int kmin, kmax;
   int projection;

   if (argc != 7) {
      printf("Usage...\n");
      printf("%s projection lon lat rad kmin kmax\n",argv[0]);
      printf("lon,lat,rad in decimal degrees\n");
      printf("projections: 0=TSC, 1=CSC, 2=QSC, 3=HSC\n");
      printf("kmin,kmax are the max and min squid resolutions\n");
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
   lon=atof(argv[2]);
   lat=atof(argv[3]);
   rad=atof(argv[4]);
   kmin=atoi(argv[5]);
   kmax=atoi(argv[6]);

   lon=lon*DD2R;
   lat=lat*DD2R;
   rad=rad*DD2R;

   if (cone_search(projection, lon, lat, rad, kmin, kmax, 
            &nfull, &full_tiles, &npart, &part_tiles) < 0) {
      fprintf(stderr,"cone_search failed in %s\n",argv[0]);
      exit(-1);
   }

   printf("%ld full_tiles:\n",nfull);
   for (i=0; i<nfull; i++) {
      printf("%ld ",(long)full_tiles[i]);
   }
   printf("\n");

   printf("%ld part_tiles:\n",npart);
   for (i=0; i<npart; i++) {
      printf("%ld ",(long)part_tiles[i]);
   }
   printf("\n");

   return(0);

}
