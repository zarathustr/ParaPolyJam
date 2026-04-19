cmake_minimum_required(VERSION 3.16)

project(test_gxy_monomial LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
set(POLYJAM_TEST_SEEDS 100 CACHE STRING "Number of synthetic RNG seeds to evaluate in the test.")

add_compile_options(-w -O3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -O3")

set(GXY_MONOMIAL_CORE_SOLVER_NAME gxy_monomial_core)
set(GXY_MONOMIAL_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../${GXY_MONOMIAL_CORE_SOLVER_NAME}")
set(GXY_MONOMIAL_CORE_SOLVER_CPP "${GXY_MONOMIAL_CORE_SOLVER_DIR}/${GXY_MONOMIAL_CORE_SOLVER_NAME}.cpp")
set(GXY_MONOMIAL_CORE_SOLVER_HPP "${GXY_MONOMIAL_CORE_SOLVER_DIR}/${GXY_MONOMIAL_CORE_SOLVER_NAME}.hpp")

find_package(Eigen3 QUIET NO_MODULE)
if(Eigen3_FOUND)
  set(POLYJAM_EIGEN_TARGET Eigen3::Eigen)
else()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")
  find_package(Eigen REQUIRED)
  add_library(polyjam_eigen INTERFACE)
  target_include_directories(polyjam_eigen INTERFACE "${EIGEN_INCLUDE_DIRS}" "${EIGEN_INCLUDE_DIR}")
  set(POLYJAM_EIGEN_TARGET polyjam_eigen)
endif()

if(NOT EXISTS "${GXY_MONOMIAL_CORE_SOLVER_CPP}" OR
   NOT EXISTS "${GXY_MONOMIAL_CORE_SOLVER_HPP}")
  message(FATAL_ERROR
    "The generated ${GXY_MONOMIAL_CORE_SOLVER_NAME}.cpp/.hpp files are missing. "
    "Build the top-level generate_gxy_monomial_solver target first.")
endif()

add_executable(test_gxy_monomial
  test_gxy_monomial.cpp
  gxy_monomial.cpp
  "${GXY_MONOMIAL_CORE_SOLVER_CPP}")

target_compile_definitions(test_gxy_monomial
  PRIVATE POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

target_include_directories(test_gxy_monomial
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${GXY_MONOMIAL_CORE_SOLVER_DIR}")

target_link_libraries(test_gxy_monomial PRIVATE ${POLYJAM_EIGEN_TARGET})
