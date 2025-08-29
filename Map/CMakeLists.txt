# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Map)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(Map Map.C)
target_link_libraries(Map Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
set(SOURCE_DIR "/home/lirui/garfield++/TempTest/GemWithRim")

foreach(_file dielectrics.dat field_triple_LH_510_corrected.txt mesh_triple_LH_510_corrected.mphtxt)
  configure_file("${SOURCE_DIR}/${_file}" ${_file} COPYONLY)
endforeach()