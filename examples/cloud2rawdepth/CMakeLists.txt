cmake_minimum_required(VERSION 3.8)
project(cloud2rawdepth)

set(CMAKE_CXX_STANDARD 17)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package (OpenCV REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/../../include)

add_executable( cloud2rawdepth 
			${PROJECT_SOURCE_DIR}/main_cloud2rawdepth.cpp)

add_library(lz4lib SHARED IMPORTED)
 #set_property(TARGET lz4lib PROPERTY IMPORTED_LOCATION "/usr/lib/x86_64-linux-gnu/liblz4.so")
 set_property(TARGET lz4lib PROPERTY IMPORTED_LOCATION "/usr/lib64/liblz4.so")
			
target_link_libraries(cloud2rawdepth 
                      Eigen3::Eigen 
                      ${OpenCV_LIBS}
                      ${PROJECT_SOURCE_DIR}/../../build/libcalibr.a
                      lz4lib
                      )
 
