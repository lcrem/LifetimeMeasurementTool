# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(lifetime)
set(libname "lifetime")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

file(GLOB HEADER_FILES "Raw*.h" "Useful*.h")

file(GLOB_RECURSE SOURCE_FILES Raw*.cxx Useful*.cxx)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
set(DICTNAME G__${libname})

ROOT_GENERATE_DICTIONARY(${DICTNAME} ${HEADER_FILES} LINKDEF LinkDef.h)

#---Create a shared library with geneated dictionary
add_library(${libname} SHARED ${SOURCE_FILES} ${DICTNAME}.cxx)
target_link_libraries(${libname} ${ROOT_LIBRARIES})


if( ${ROOT_VERSION} VERSION_GREATER "5.99.99")
    add_custom_target(${DICTNAME}.pcm DEPENDS ${DICTNAME})
endif()

message("DICTNAME is set to ${DICTNAME}" )
message("HEADER_FILES is set to ${HEADER_FILES}" )
message("SOURCE_FILES is set to ${SOURCE_FILES}" )

#---Create  a main program using the library
add_executable(makeTGraphs makeTGraphs.cxx) 
target_link_libraries(makeTGraphs ${ROOT_LIBRARIES} lifetime)

add_executable(findAllAverages findAllAverages.cxx) 
target_link_libraries(findAllAverages ${ROOT_LIBRARIES} lifetime)

add_executable(calculateLifetimeWithScope calculateLifetimeWithScope.cxx) 
target_link_libraries(calculateLifetimeWithScope ${ROOT_LIBRARIES} lifetime)



################### DOXYGEN STUFF

find_package(Doxygen)
option(BUILD_DOCUMENTATION "Create and install the HTML based API        
documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    if(NOT DOXYGEN_FOUND)
         message(FATAL_ERROR "Doxygen is needed to build the documentation.")
    endif()

    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen.config)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    message("Doxygen build started.")

    add_custom_target(doc
                      COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile_in}
                      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                      COMMENT "Generating API documentation with Doxygen"
                      VERBATIM)

    #    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION     share/doc)
endif()
