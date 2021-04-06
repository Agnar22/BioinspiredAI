cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER "/usr/bin/gcc-9")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-9")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH "/usr/include/opencv/build")

project(MOEA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -pthread -lpthread -O3 -march=native")

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_library(ga dir.h prim.h prim.cpp file.h file.cpp individual.h individual.cpp objective.h)
add_executable(main main.cpp)
target_link_libraries(main ga ${OpenCV_LIBS})