cmake_minimum_required(VERSION 2.8.12)
project(mpm_gmsh_generator)

#compiler flags (through conversation with others, not sure what most of them do.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -O3 -fPIC -g -D_GLIBCXX_USE_CXX11_ABI=0")

#set main source file
set(SOURCE_FILES main.cpp)

add_executable(mpm_gmsh_generator ${SOURCE_FILES})
