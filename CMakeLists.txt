cmake_minimum_required(VERSION 2.8.12)
project(mpm_3d)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/driver/main.cpp)
set(SOURCE_FILES ${SOURCE_FILES}
        src/objects/body.cpp
        src/functions/process.cpp
        src/functions/spmd.cpp)

include_directories(include)

add_executable(mpm_3d ${SOURCE_FILES})