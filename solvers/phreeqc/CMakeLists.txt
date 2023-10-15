cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(runFreak)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#for torch lpthread sugges tot add -DDEAL_II_LINKER_FLAGS=-lpthread, does not work
set(DEFINITIONS_COMPILE "-std=c++14 -D_GLIBCXX_USE_CXX11_ABI=0 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wno-unused-parameter -Wno-overloaded-virtual -Wno-unused-variable -Wno-unused-local-typedef -Wno-invalid-offsetof -Wno-deprecated-register -Wno-undefined-var-template -O0 -g -DFULLDEBUG -DNoRepository -ftemplate-depth-100 -fPIC")
# =====================================================================================================
# Compiling configure
add_definitions("${DEFINITIONS_COMPILE}")
add_definitions("-Dlinux64")
# ==========================================

include_directories(/usr/local/lib) # for phreeqc
include_directories(/usr/lib) # for torch libs

#link_directories(${FOAM_LIB_DIR} ${FOAM_USER_LIBBIN} ${FOAM_LIB_DIR} ${PATH_LIB_OPENMPI})
link_directories(/home/oatteia/lib)

add_executable(runFreak runFreak.cpp)
#target_link_libraries(muFlowRTnn "${MPI_C_LIBRARIES}" "${TORCH_LIBRARIES}")
target_link_libraries(${PROJECT_NAME} PhreeqcRM)