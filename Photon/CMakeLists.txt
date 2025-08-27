# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Photon)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(Photon Photon.C)
target_link_libraries(Photon Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
set(SOURCE_DIR "/home/.../GemWithRim")

foreach(_file dielectrics.dat field_single_gem_560_190.txt mesh_single_gem_560_190.mphtxt)
  configure_file("${SOURCE_DIR}/${_file}" ${_file} COPYONLY)
endforeach()