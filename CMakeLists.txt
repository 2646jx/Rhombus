cmake_minimum_required(VERSION 3.13)
project(Rhombus LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS Off)

add_compile_options(-fPIC)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message(STATUS "CMAKE_BUILD_TYPE not detected, set to default (${CMAKE_BUILD_TYPE})")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(RHOMBUS_DEBUG ON)
else()
    set(RHOMBUS_DEBUG OFF)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(RHOMBUS_BUILD_TESTS "Enable tests" OFF)

include(CMakeFindDependencyMacro)

set(CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
find_package(SEAL REQUIRED)
if(SEAL_FOUND)
    message(STATUS "SEAL found!")
else()
    message(FATAL_ERROR "SEAL not found")
endif()

message(STATUS "SEAL_INCLUDE_DIR" ${SEAL_INCLUDE_DIR})

add_subdirectory(src)

if(RHOMBUS_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
