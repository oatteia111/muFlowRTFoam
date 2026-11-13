cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(muFlowRT)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION")
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/../VERSION" MUFLOWRT_VERSION)
string(STRIP "${MUFLOWRT_VERSION}" MUFLOWRT_VERSION)
if(MUFLOWRT_VERSION STREQUAL "")
    message(FATAL_ERROR "VERSION must not be empty")
endif()
set(MUFLOWRT_RELEASE_TAG "$ENV{MUFLOWRT_RELEASE_TAG}")
if(MUFLOWRT_RELEASE_TAG AND NOT MUFLOWRT_RELEASE_TAG STREQUAL "v${MUFLOWRT_VERSION}")
    message(FATAL_ERROR "Release tag ${MUFLOWRT_RELEASE_TAG} does not match VERSION (${MUFLOWRT_VERSION})")
endif()
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/version.hpp" @ONLY)

find_package( MPI REQUIRED )
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
# add this line only when you are using openmpi which has a different c++ bindings

#list(APPEND CMAKE_PREFIX_PATH "/home/oatteia/muFlowRT/pytorch1/libtorch")
#find_package(Torch REQUIRED)

# =====================================================================================
#             OpenFOAM configurations
# -------------------------------------------------------------------------------------
#Check valid OpenFOAM
if(DEFINED ENV{WM_PROJECT_DIR})
	MESSAGE(STATUS "OpenFOAM: " $ENV{WM_PROJECT_DIR})
else()
	message(FATAL_ERROR "The OpenFOAM bashrc is not sourced")
endif(DEFINED ENV{WM_PROJECT_DIR})
set(FOAM_VERSION $ENV{WM_PROJECT_VERSION}) 
set(FOAM_DIR $ENV{WM_PROJECT_DIR})
set(FOAM_LIB_DIR $ENV{FOAM_LIBBIN})
set(FOAM_USER_LIBBIN $ENV{FOAM_USER_LIBBIN})
set(FOAM_SRC $ENV{FOAM_SRC})
# lines set and find from https://github.com/stuart23/cmake-OpenFOAM/blob/master/CMakeLists.txt
set( CMAKE_MODULE_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )
set( CMAKE_SHARED_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )
set( CMAKE_EXE_LINKER_FLAGS "-Xlinker --copy-dt-needed-entries -Xlinker --no-as-needed" )
find_package( MPI REQUIRED )
#find_package( SCOTCH REQUIRED )
#find_package( PTSCOTCH REQUIRED )
find_package( FLEX REQUIRED )
find_package( ZLIB REQUIRED )
find_package( Boost COMPONENTS thread REQUIRED )
#find_package( METIS REQUIRED )
#find_package( Zoltan REQUIRED )
find_package( Threads REQUIRED )
find_package( CGAL REQUIRED )


#set(PATH_LIB_OPENMPI "openmpi-system")  # Foundation version #to put below in case -D_GLIBCXX_USE_CXX17_ABI=0
#for torch lpthread sugges tot add -DDEAL_II_LINKER_FLAGS=-lpthread, does not work
set(DEFINITIONS_COMPILE "-std=c++11 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedef -Wno-invalid-offsetof -Wno-deprecated-register -Wno-undefined-var-template -O0 -g -DFULLDEBUG -DNoRepository -ftemplate-depth-100 -fPIC")
# =====================================================================================================
# Compiling configure
add_definitions("${DEFINITIONS_COMPILE}")
add_definitions(-DUSE_MPI)
add_definitions("-Dlinux64")
# ==========================================
include_directories(. 
                    ${FOAM_SRC}/engine/lnInclude  
                    ${FOAM_SRC}/ODE/lnInclude  
                    ${FOAM_SRC}/OpenFOAM/lnInclude  
					${FOAM_SRC}/OSspecific/POSIX/lnInclude  
                    ${FOAM_SRC}/fvOptions/lnInclude 
					${FOAM_SRC}/finiteVolume/lnInclude 
                    ${FOAM_SRC}/meshTools/lnInclude
					) #

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libraries/phaseModels/lnInclude)
include_directories(usr/local/lib/include) # for phreeqc
include_directories(/usr/local/lib) # for phreeqc
include_directories(/usr/lib) # for torch libs
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(${MPI_CXX_INCLUDE_PATH})

set(LINK_LIBS ${LINK_LIBS} ${MPI_CXX_LIBRARIES})
link_directories(${FOAM_LIB_DIR} ${FOAM_USER_LIBBIN} ${FOAM_LIB_DIR} ${PATH_LIB_OPENMPI})

add_executable(${PROJECT_NAME} muFlowRT.C)
target_include_directories(${PROJECT_NAME} BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(${PROJECT_NAME} PhreeqcRM ${MPI_CXX_LIBRARIES})
#target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} pthread OpenFOAM dl m finiteVolume meshTools phase)
