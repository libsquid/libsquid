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

project(LIBSQUID)
cmake_minimum_required(VERSION 2.8)

# Read libsquid.h and get version and release string
file(READ "${LIBSQUID_SOURCE_DIR}/libsquid.h" LIBSQUID_H)
string(FIND "${LIBSQUID_H}" "LIBSQUID_VERSION" LIBSQUID_VERSION_CHECK)
if (LIBSQUID_VERSION_CHECK EQUAL -1)
   message(FATAL_ERROR "version number not in libsquid.h")
endif()
string(REGEX REPLACE ".*#define LIBSQUID_VERSION[^0-9]*([0-9]+)\\.([0-9]+)\\.([0-9]+).*" 
	     "\\1.\\2.\\3" LIBSQUID_VERSION_STRING "${LIBSQUID_H}")
if (LIBSQUID_VERSION_STRING STREQUAL ifdef)
   message(FATAL_ERROR "version string not set in libsquid.h")
endif()
string(FIND "${LIBSQUID_H}" "LIBSQUID_RELEASE" LIBSQUID_RELEASE_CHECK)
if (LIBSQUID_RELEASE_CHECK EQUAL -1)
   set(LIBSQUID_RELEASE_STRING "")
else()
   string(REGEX REPLACE ".*#define LIBSQUID_RELEASE[^0-9a-zA-Z\\.]*([0-9a-zA-Z\\.\\-]+).*" 
   	     "\\1" LIBSQUID_RELEASE_STRING "${LIBSQUID_H}")
   message("LIBSQUID_VERSION_STRING ${LIBSQUID_VERSION_STRING}")
   if (LIBSQUID_RELEASE_STRING STREQUAL ifdef)
      set(LIBSQUID_RELEASE_STRING "")
   endif()
   message("LIBSQUID_RELEASE_STRING ${LIBSQUID_RELEASE_STRING}")
endif()

# Get list of source and header files
set(CMAKE_INCLUDE_CURRENT_DIR ON)
file(GLOB LIBSQUID_SOURCES *.c)
file(GLOB LIBSQUID_HEADERS *.h)

# Add static and shared library targets
add_library(libsquid_static STATIC ${LIBSQUID_SOURCES})
set_target_properties(libsquid_static PROPERTIES OUTPUT_NAME squid)
add_library(libsquid_shared SHARED ${LIBSQUID_SOURCES})
set_target_properties(libsquid_shared PROPERTIES 
			     OUTPUT_NAME squid
			     SOVERSION ${LIBSQUID_VERSION_STRING})
add_subdirectory(bin)

install(TARGETS libsquid_shared
	LIBRARY	DESTINATION lib
	ARCHIVE	DESTINATION lib
	COMPONENT libsquid
)
install(FILES ${LIBSQUID_HEADERS}
	DESTINATION include
	COMPONENT libsquid
)
install(TARGETS libsquid_static
	LIBRARY	DESTINATION lib
	ARCHIVE	DESTINATION lib
	COMPONENT libsquid
)

add_custom_target(libsquid_dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source )

set(CPACK_PACKAGE_NAME "libsquid")
set(CPACK_PACKAGE_VERSION ${LIBSQUID_VERSION_STRING})
if (LIBSQUID_RELEASE_STRING STREQUAL "")
   set(CPACK_PACKAGE_RELEASE "")
else()
   set(CPACK_PACKAGE_RELEASE "~${LIBSQUID_RELEASE_STRING}")
endif()
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
    "The Spherical-Cube Quad-Tree Unique Index Library (LibSQUID)")
set(CPACK_PACKAGE_CONTACT "Jim Wren; jwren@lanl.gov")
set(CPACK_PACKAGE_VENDOR ${CPACK_PACKAGE_CONTACT})
set(CPACK_RPM_PACKAGE_RELEASE "${LIBSQUID_RELEASE_STRING}")
set(CPACK_GENERATOR "RPM" "DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "~$;[.]bak$/;[.]swp$;[.]svn$")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/packages)
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}${CPACK_PACKAGE_RELEASE}.${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}${CPACK_PACKAGE_RELEASE}")
include(CPack)
