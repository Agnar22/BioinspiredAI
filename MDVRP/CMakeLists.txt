cmake_minimum_required(VERSION 3.10)


set(CMAKE_C_COMPILER "/usr/bin/gcc-9")
set(CMAKE_CXX_COMPILER "/usr/bin/g++-9")
set(CMAKE_CXX_STANDARD 17)

project(MDVRP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -pthread -lpthread -O3 -march=native")


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_executable(main main.cpp ga.h ga.cpp individual.h individual.cpp problem.h problem.cpp file.h file.cpp)