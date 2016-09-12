cmake_minimum_required(VERSION 2.8.12)
project(mpm_3d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/driver/main.cpp)
set(SOURCE_FILES ${SOURCE_FILES}
        src/objects/body.cpp
        src/objects/material1.cpp
        src/objects/material2.cpp
        src/functions/process.cpp
        src/functions/tensor.cpp)

include_directories(include)
include_directories(src/materials)

add_executable(mpm_3d ${SOURCE_FILES})