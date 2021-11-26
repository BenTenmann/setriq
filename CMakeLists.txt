cmake_minimum_required(VERSION 3.17)
set(PROJECT_NAME setriq)
set(PROJECT_SOURCE_DIR src)
set(PROJECT_TEST_DIR tests)

file(
        GLOB_RECURSE
        SOURCE
        ${PROJECT_SOURCE_DIR}/*.cpp
        ${PROJECT_SOURCE_DIR}/*.h
)

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 14)
include_directories(include/setriq lib/jsoncpp)

add_executable(${PROJECT_NAME} ${SOURCE} ${PROJECT_TEST_DIR}/main.cpp)

find_package(pybind11 CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE mylib pybind11::lto pybind11::embed pybind11::module)
