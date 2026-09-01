cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(runFreak2)
set(DEFINITIONS_COMPILE "-std=c++11 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedef -Wno-invalid-offsetof -Wno-deprecated-register -Wno-undefined-var-template -O0 -g -DFULLDEBUG -DNoRepository -ftemplate-depth-100 -fPIC")
# =====================================================================================================
# Compiling configure
add_definitions("${DEFINITIONS_COMPILE}")
add_definitions("-Dlinux64")
# ==========================================
include_directories(/usr/local/include) # for phreeqc
#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(${PROJECT_NAME} runFreak2.cpp)
target_link_libraries(${PROJECT_NAME} dl m phreeqcrm)