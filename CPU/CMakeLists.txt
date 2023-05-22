cmake_minimum_required(VERSION 3.14)
project(CPU)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-g")

include_directories(.)
include_directories(Common)
include_directories(Multiple)
include_directories(Struct)

add_executable(CPU main.cpp)
