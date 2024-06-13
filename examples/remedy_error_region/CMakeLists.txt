cmake_minimum_required(VERSION 3.8)
project(remedy_error_region)

set(CMAKE_CXX_STANDARD 17)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package (OpenCV REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/../../include)

add_executable(remedy_error_region 
			${PROJECT_SOURCE_DIR}/main_remedy_error_region.cpp)

add_library(flannlib SHARED IMPORTED)
set_property(TARGET flannlib PROPERTY IMPORTED_LOCATION "/usr/lib64/libflann_cpp.so")

add_library(lz4lib SHARED IMPORTED)
set_property(TARGET lz4lib PROPERTY IMPORTED_LOCATION "/usr/lib64/liblz4.so")
			
target_link_libraries(remedy_error_region  
                      Eigen3::Eigen 
                      ${OpenCV_LIBS}
                      ${PROJECT_SOURCE_DIR}/../../build/libcalibr.a
                      flannlib
                      lz4lib
                      )
 
