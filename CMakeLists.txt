cmake_minimum_required(VERSION 3.17)
set(PROJECT_NAME setriq)
set(PROJECT_SOURCE_DIR setriq/_C)
set(PROJECT_TEST_DIR tests)

file(
        GLOB_RECURSE
        SOURCE
        "${PROJECT_SOURCE_DIR}/*.cpp"
        "${PROJECT_SOURCE_DIR}/*.h"
)

get_filename_component(PYTHON_WRAPPER "${PROJECT_SOURCE_DIR}/pythonWrapper.cpp" ABSOLUTE)
list(REMOVE_ITEM SOURCE "${PYTHON_WRAPPER}")

project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 14)
include_directories(include/setriq)

add_executable(${PROJECT_NAME} ${SOURCE} ${PROJECT_TEST_DIR}/main.cpp)

find_package(jsoncpp CONFIG REQUIRED)
target_link_libraries(${PROJECT_NAME} PRIVATE jsoncpp_object jsoncpp_static JsonCpp::JsonCpp)