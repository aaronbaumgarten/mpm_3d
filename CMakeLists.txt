cmake_minimum_required(VERSION 2.8.12)
project(mpm_v2)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -O3 -fPIC -D_GLIBCXX_USE_CXX11_ABI=0")

# dlopen required for file read
find_package(DL REQUIRED)

set(SOURCE_FILES src/main/main.cpp)
set(SOURCE_FILES ${SOURCE_FILES}
        src/objects/job.cpp
        src/objects/serializer.cpp
        src/objects/driver.cpp
        src/objects/solver.cpp
        src/objects/body.cpp
        src/objects/contact.cpp
        src/objects/grid.cpp
        src/objects/config.cpp
        src/objects/points.cpp
        src/objects/nodes.cpp
        src/objects/boundary.cpp
        src/objects/material.cpp
        src/functions/stringparser.cpp
        )

include_directories(include)
include_directories(/usr/include/eigen3)
include_directories(${DL_INCLUDE_DIR})

add_subdirectory(src/serializers)
add_subdirectory(src/drivers)
add_subdirectory(src/solvers)
add_subdirectory(src/grids)
add_subdirectory(src/materials)
add_subdirectory(src/boundaries)
add_subdirectory(src/contacts)
add_executable(mpm_v2 ${SOURCE_FILES})

target_link_libraries(mpm_v2 ${DL_LIBRARIES})

#set(UNIT_TEST_SOURCE_FILES
#    src/tests/unitTestMain.cpp
#    src/tests/unitShapeFunctions.cpp)
#add_executable(unit ${UNIT_TEST_SOURCE_FILES})
