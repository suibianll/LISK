# Project properties
project(lisk)

cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 20)

# #  -fsanitize=address
set(CMAKE_CXX_FLAGS "-std=c++17 -faligned-new -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

find_library(TBB_LIBRARIES NAMES tbb)
find_package(OpenMP)
# add_compile_options(-faligned-new -march=native -g)
include_directories(src)

# Now add what we want to compile
add_subdirectory(src)

add_executable(microbench test.cpp)
target_link_libraries(microbench lisk glog jemalloc)

add_executable(loadbench testload.cpp)
target_link_libraries(loadbench lisk glog jemalloc)