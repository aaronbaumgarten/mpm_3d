cmake_minimum_required(VERSION 2.8.12)
project(mpm_3d)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -O3")

# dlopen required for file read
find_package(DL REQUIRED)

set(SOURCE_FILES src/driver/main.cpp)
set(SOURCE_FILES ${SOURCE_FILES}
        src/objects/body.cpp
        src/objects/material.cpp
        src/objects/boundary.cpp
        src/objects/contact.cpp
        src/objects/node.cpp
        src/objects/particle.cpp
        src/objects/element.cpp
        src/objects/mpmio.cpp
        src/functions/process.cpp
        src/functions/tensor.cpp
        src/functions/loading.cpp
        src/tests/test.cpp)

include_directories(include)
include_directories(/usr/include/eigen3)
include_directories(${DL_INCLUDE_DIR})

add_subdirectory(src/materials)
add_subdirectory(src/boundaries)
add_subdirectory(src/contacts)
add_executable(mpm_3d ${SOURCE_FILES})

target_link_libraries(mpm_3d ${DL_LIBRARIES})

#set(UNIT_TEST_SOURCE_FILES
#    src/tests/unitTestMain.cpp
#    src/tests/unitShapeFunctions.cpp)
#add_executable(unit ${UNIT_TEST_SOURCE_FILES})
