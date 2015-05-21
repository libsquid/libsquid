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
   int k, face, doprint=0;
   long xl,yl,nside;
   double xd,yd;
   squid_type squid;
   double lon1,lat1,lon2,lat2;
   int face1,face2;
   double x1,y1,x2,y2;
   int projection;
   double tol; // tolerance between floats

   if (argc != 2) {
      printf("Example usage...\n");
      printf("%s projection\n",argv[0]);
      printf("projections: 0=TSC, 1=CSC, 2=QSC, 3=HSC\n");
      exit(-1);
   }
   projection=atoi(argv[1]);
   tol=0.00001;
   if (projection == TSC) {
      printf("TSC Projection\n");
   } else if (projection == CSC) {
      printf("CSC Projection\n");
      tol=0.001;
   } else if (projection == QSC) {
      printf("QSC Projection\n");
   } else if (projection == HSC) {
      printf("HSC Projection\n");
   } else {
      printf("Unknown projection! Using HSC.\n");
      projection=HSC;
   }

   k=10;
   nside=(long)pow(2,k);
   printf("nside=%ld\n",nside);

   for (face=0; face<6; face++) {
      printf("face=%d\n",face);
      for (yl=0; yl<nside; yl++) {
         for (xl=0; xl<nside; xl++) {
            xd=((double)xl+0.5)/(double)nside;
            yd=((double)yl+0.5)/(double)nside;
            xyf2sph(projection,xd,yd,face,&lon1,&lat1);
            sph2xyf(projection,lon1,lat1,&x1,&y1,&face1);
            sph2squid(projection,lon1,lat1,k,&squid);
            squid2sph(projection,squid,&lon2,&lat2);
            sph2xyf(projection,lon2,lat2,&x2,&y2,&face2);

            if (fabs(lon1-lon2) > tol) {
               printf("lon1 != lon2\n");
               doprint=1;
            }
            if (fabs(lat1-lat2) > tol) {
               printf("lat1 != lat2\n");
               doprint=1;
            }
            if (fabs(xd-x1) > tol) {
               printf("xd != x1\n");
               doprint=1;
            }
            if (fabs(yd-y1) > tol) {
               printf("yd != y1\n");
               doprint=1;
            }
            if (face != face1) {
               printf("facel != face1\n");
               doprint=1;
            }
            if (fabs(x2-x1) > tol) {
               printf("x2 != x1\n");
               doprint=1;
            }
            if (fabs(y2-y1) > tol) {
               printf("y2 != y1\n");
               doprint=1;
            }
            if (face2 != face1) {
               printf("face2 != face1\n");
               doprint=1;
            }
            if (doprint == 1) break;
         }
         if (doprint == 1) break;
      }
      if (doprint == 1) break;
   }

   if (doprint) {
      printf("xl=%ld yl=%ld\n",xl,yl);
      printf("xd=%.10f yd=%.10f face=%d\n",xd,yd,face);
      printf("lon1=%.10f lat1=%.10f\n",lon1/DD2R,lat1/DD2R);
      printf("x1=%.3f y1=%.3f face1=%d\n",x1,y1,face1);
      printf("squid=%ld\n",(long)squid);
      printf("lon2=%.5f lat2=%.5f\n",lon2/DD2R,lat2/DD2R);
      printf("x2=%.3f y2=%.3f face2=%d\n",x2,y2,face2);
   }

   return(0);

}


