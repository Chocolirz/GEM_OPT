# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(Field)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(Field Field.C)
target_link_libraries(Field Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
set(SOURCE_DIR "/home/.../GemWithRim")

foreach(_file dielectrics.dat field_single_gem_580_170.txt mesh_single_gem_580_170.mphtxt)
  configure_file("${SOURCE_DIR}/${_file}" ${_file} COPYONLY)
endforeach()