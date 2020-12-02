# CMakeLists for RADGEN
project(RADGEN)

# Set up
cmake_minimum_required (VERSION 3.0)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

set(CMAKE_CXX_FLAGS "-std=c++11")

enable_language(Fortran)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)

include_directories(${INCLUDES})

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.cc ${CMAKE_SOURCE_DIR}/src/fortran/*.f)
file(GLOB HEADERS ${CMAKE_SOURCE_DIR}/include/*.hh)

add_executable(run_radgen run_radgen.cc ${SOURCES} ${HEADERS})
target_link_libraries(run_radgen ${ROOT_LIBRARIES})
