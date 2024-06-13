cmake_minimum_required(VERSION 3.8)
project(calibr)

set(CMAKE_CXX_STANDARD 17)

find_package (Eigen3 REQUIRED NO_MODULE)
find_package (OpenCV REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include
	/usr/local/include/opencv4)

add_library(calibr 
			${PROJECT_SOURCE_DIR}/src/calibr.cpp
			${PROJECT_SOURCE_DIR}/src/homography2d.cpp
			${PROJECT_SOURCE_DIR}/src/image.cpp
			${PROJECT_SOURCE_DIR}/src/numer_math.cpp
			${PROJECT_SOURCE_DIR}/src/reconst.cpp
            ${PROJECT_SOURCE_DIR}/src/cloud2.cpp
			${PROJECT_SOURCE_DIR}/src/knn_search.cpp
			${PROJECT_SOURCE_DIR}/src/search_grid.cpp)

#add_library(lz4lib SHARED IMPORTED)
# set_property(TARGET lz4lib PROPERTY IMPORTED_LOCATION "/usr/lib64/liblz4.so")
	
target_link_libraries(calibr Eigen3::Eigen)

add_subdirectory(examples)

