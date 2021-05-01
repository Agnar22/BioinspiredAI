cmake_minimum_required(VERSION 3.10)

set(CMAKE_C_COMPILER "/usr/bin/gcc-9")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-9")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH "/usr/include/opencv/build")

project(MOEA)

enable_testing()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -pthread -lpthread -O3 -march=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -pthread -lpthread")


find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_library(individual individual.h individual.cpp objective.h objective.cpp prim.h prim.cpp dir.h)
add_library(ga ga.h ga.cpp nsga.h nsga.cpp)
add_library(file file.h file.cpp)
add_library(slic slic.h slic.cpp)
add_library(kruskal kruskal.h kruskal.cpp)

#add_library(ga ga.h ga.cpp dir.h prim.h prim.cpp nsga.h nsga.cpp file.h file.cpp individual.h individual.cpp objective.h objective.cpp)
add_executable(main main.cpp)
target_link_libraries(main ga individual file slic kruskal ${OpenCV_LIBS})

add_subdirectory(tests)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()
