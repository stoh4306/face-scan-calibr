cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(marching_cubes)

find_package(CUDA REQUIRED)
message(${CUDA_TOOLKIT_ROOT_DIR})

set(CMAKE_CUDA_COMPILER ${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc)

enable_language("CUDA")

if(NOT DEFINED CMAKE_CXX_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED True)
endif()


include_directories(${PROJECT_SOURCE_DIR})
link_directories(/usr/lib/x86_64-linux-gnu
                ${CUDA_TOOLKIT_ROOT_DIR}/lib64)

add_executable(marching_cubes
               main.cpp
               marching_cubes.cpp
               marching_cubes_cuda.cu)

target_link_libraries(marching_cubes cudart ${PROJECT_SOURCE_DIR}/../../build/libcalibr.a)