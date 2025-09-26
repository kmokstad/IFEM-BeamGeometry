cmake_minimum_required(VERSION 3.22)

option(IFEM_USE_FFALIB "Compile with FFaLib?" OFF)

project(BeamGeometry LANGUAGES C CXX)

message(STATUS "Configuring project ${PROJECT_NAME}")

if(NOT TARGET IFEM)
  if(IFEM_AS_SUBMODULE)
    add_subdirectory(../../.. IFEM)
    set(IFEM_PATH ${PROJECT_SOURCE_DIR}/../../..)
  else()
    get_filename_component(BUILD_DIR ${PROJECT_BINARY_DIR} NAME)
    set(IFEM_DIR ${PROJECT_SOURCE_DIR}/../../../${BUILD_DIR})
    find_package(IFEM REQUIRED)
  endif()
endif()

if(NOT TARGET Elasticity)
  add_subdirectory(.. Elasticity)
endif()

if(NOT TARGET Beam)
  add_subdirectory(../Beam Beam)
endif()

if(IFEM_USE_FFALIB)
  message(STATUS "Compiling with use of FFaLib")
  add_subdirectory(FFaLib)
endif()

ifem_add_library(
  NAME
    BeamGeometry
  SOURCES
    readTesselation.C
    SIMBeamGeometry.C
  HEADERS
    readTesselation.h
    SIMBeamGeometry.h
  LIBRARIES
    Beam
)
target_compile_definitions(BeamGeometry INTERFACE HAS_GEOMETRY=1)

if(TARGET FFaLib)
  target_link_libraries(BeamGeometry PUBLIC FFaLib)
endif()

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
  message(STATUS "Adding dummy testapps target for project ${CMAKE_PROJECT_NAME}")
  add_custom_target(testapps COMMAND echo "No apps here" DEPENDS BeamGeometry)
  install(CODE "message(\"No installation here\")")
endif()
