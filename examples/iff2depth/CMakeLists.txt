cmake_minimum_required(VERSION 3.8)
project(iff2depth)

set(CMAKE_CXX_STANDARD 17)

find_package (OpenCV REQUIRED)

#include_directories(${PROJECT_SOURCE_DIR}/)

add_executable(iff2depth
			${PROJECT_SOURCE_DIR}/main_iff2depth.cpp
			${PROJECT_SOURCE_DIR}/renew_iff.cpp
			${PROJECT_SOURCE_DIR}/renew_ifflib.cpp)
			
target_link_libraries(iff2depth  
                      ${OpenCV_LIBS})
 
