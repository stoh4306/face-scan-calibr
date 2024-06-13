cmake_minimum_required(VERSION 3.8)
project(polaRawToRGB)

set(CMAKE_CXX_STANDARD 17)

find_package (OpenCV REQUIRED)

#include_directories(${PROJECT_SOURCE_DIR}/)

add_executable(polaRawToRGB
			${PROJECT_SOURCE_DIR}/main.cpp)
			
target_link_libraries(polaRawToRGB  
                      ${OpenCV_LIBS})
 
