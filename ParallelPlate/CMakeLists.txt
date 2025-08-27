cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(ParallelPlate)
if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()

add_executable(ParallelPlate ParallelPlate.C)
target_link_libraries(ParallelPlate Garfield::Garfield)
