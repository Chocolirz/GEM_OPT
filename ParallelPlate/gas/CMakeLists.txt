cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(GasFile)
if(NOT TARGET Garfield::Garfield)
    find_package(Garfield REQUIRED)
endif()

add_executable(GasFile GasFile.C)
target_link_libraries(GasFile Garfield::Garfield)
