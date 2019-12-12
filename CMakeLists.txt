cmake_minimum_required(VERSION 2.8.12)
project(mpm_v3)

#compiler flags (through conversation with others, not sure what most of them do.
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -O3 -fPIC -g -D_GLIBCXX_USE_CXX11_ABI=0 -pthread")

#set main source file
set(SOURCE_FILES src/main/main.cpp)

#set include directories
include_directories(include)
include_directories(src)
include_directories(/usr/include/)
include_directories(/usr/include/eigen3)
include_directories(/usr/local/include/)
include_directories(/usr/local/include/eigen3) #eigen is required for this

add_subdirectory(src/tests)
add_subdirectory(src/registry)
add_subdirectory(src/algebra)

add_subdirectory(src/config)

add_subdirectory(src/fvm)

add_subdirectory(src/objects/serializers)
add_subdirectory(src/objects/drivers)
add_subdirectory(src/objects/solvers)

add_subdirectory(src/objects/bodies)
add_subdirectory(src/objects/contacts)
add_subdirectory(src/objects/grids)

add_subdirectory(src/objects/points)
add_subdirectory(src/objects/nodes)
add_subdirectory(src/objects/materials)
add_subdirectory(src/objects/boundaries)

add_executable(mpm_v3 ${SOURCE_FILES})
