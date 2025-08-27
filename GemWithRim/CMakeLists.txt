# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(GemWithRim)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(GemWithRim GemWithRim.C)
target_link_libraries(GemWithRim Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file dielectrics.dat field_single_gem_600_210.txt mesh_single_gem_600_210.mphtxt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()