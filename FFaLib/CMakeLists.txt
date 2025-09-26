cmake_minimum_required(VERSION 3.22)

project(FFaLib LANGUAGES C CXX)

ifem_add_library(
  NAME
    FFaLib
  SOURCES
    FFaAlgebra/FFaBody.C
    FFaAlgebra/FFaBodyParser.C
    FFaAlgebra/FFaMat33.C
    FFaAlgebra/FFaMat34.C
    FFaAlgebra/FFaMath.C
    FFaAlgebra/FFaTensor1.C
    FFaAlgebra/FFaTensor2.C
    FFaAlgebra/FFaTensor3.C
    FFaAlgebra/FFaTensorTransforms.C
    FFaAlgebra/FFaVec3.C
  HEADERS
    FFaAlgebra/FFaBody.H
    FFaAlgebra/FFaMat33.H
    FFaAlgebra/FFaMat34.H
    FFaAlgebra/FFaMath.H
    FFaAlgebra/FFaTensor1.H
    FFaAlgebra/FFaTensor2.H
    FFaAlgebra/FFaTensor3.H
    FFaAlgebra/FFaTensorTransforms.H
    FFaAlgebra/FFaVec3.H
)

target_compile_definitions(FFaLib INTERFACE HAS_FFALIB=1)
target_include_directories(FFaLib PRIVATE ${PROJECT_SOURCE_DIR}/..)
