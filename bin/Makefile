# -------------------------- LICENSE -----------------------------------
#
#  This file is part of the LibSQUID software libraray.
#
#  LibSQUID is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published
#  by the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  LibSQUID is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright 2014 James Wren and Los Alamos National Laboratory
#

TARGET_BINS = test_libsquid squidinfo sph2squid tile_sph2xy tile_xy2sph \
	cone_search
GCC     = gcc
CFLAGS  = -g -Wall -fPIC -I../
LDFLAGS = -L../ -Wl,-Bstatic -lsquid -Wl,-Bdynamic -lm

all: $(TARGET_BINS)

%.o: %.c
	$(GCC) -c $(CFLAGS) -o $@ $<

$(TARGET_BINS): % : %.o
	$(GCC) $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -f *.o $(TARGET_BINS)

