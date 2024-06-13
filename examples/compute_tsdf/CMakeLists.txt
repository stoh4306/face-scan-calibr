cmake_minimum_required(VERSION 3.8)
project(compute_tsdf)

set(CMAKE_CXX_STANDARD 17)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package (OpenCV REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/../../include)

add_executable( compute_tsdf 
			${PROJECT_SOURCE_DIR}/main_compute_tsdf.cpp)

add_library(lz4lib SHARED IMPORTED)
 set_property(TARGET lz4lib PROPERTY IMPORTED_LOCATION "/usr/lib64/liblz4.so")
			
target_link_libraries(compute_tsdf  
                      Eigen3::Eigen 
                      ${OpenCV_LIBS}
                      ${PROJECT_SOURCE_DIR}/../../build/libcalibr.a
                      lz4lib
)
 
