# -------------------------- LICENSE -----------------------------------
#
# This file is part of the LibSQUID software library.
#
# LibSQUID is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LibSQUID is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with LibSQUID.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2014 James Wren and Los Alamos National Laboratory
#

cmake_minimum_required(VERSION 2.8)
file(GLOB LIBSQUID_EXE RELATIVE ${LIBSQUID_SOURCE_DIR}/bin *.c)
include_directories(${LIBSQUID_SOURCE_DIR})
link_directories(${LIBSQUID_BINARY_DIR})
link_libraries(squid m)

foreach(infile ${LIBSQUID_EXE})
   string(REPLACE ".c" "" outfile ${infile})
   add_executable(${outfile} ${infile})
endforeach(infile)
