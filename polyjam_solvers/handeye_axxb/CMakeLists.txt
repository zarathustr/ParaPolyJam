cmake_minimum_required(VERSION 3.16)

project(test_handeye_axxb LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

add_compile_options(-w -O3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -O3")

set(POLYJAM_HAND_EYE_PAIRS 2 CACHE STRING "Number of A/B motion pairs expected by the generated solver.")
set(POLYJAM_TEST_SEEDS 100 CACHE STRING "Number of synthetic RNG seeds to evaluate in the test.")
set(HAND_EYE_CORE_SOLVER_NAME handeye_axxb_q)
set(HAND_EYE_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../${HAND_EYE_CORE_SOLVER_NAME}")
set(HAND_EYE_CORE_SOLVER_CPP "${HAND_EYE_CORE_SOLVER_DIR}/${HAND_EYE_CORE_SOLVER_NAME}.cpp")
set(HAND_EYE_CORE_SOLVER_HPP "${HAND_EYE_CORE_SOLVER_DIR}/${HAND_EYE_CORE_SOLVER_NAME}.hpp")

find_package(Eigen3 QUIET NO_MODULE)
if(Eigen3_FOUND)
  set(POLYJAM_EIGEN_TARGET Eigen3::Eigen)
else()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
  find_package(Eigen REQUIRED)
  add_library(polyjam_eigen INTERFACE)
  target_include_directories(polyjam_eigen INTERFACE "${EIGEN_INCLUDE_DIRS}" "${EIGEN_INCLUDE_DIR}")
  set(POLYJAM_EIGEN_TARGET polyjam_eigen)
endif()

if(NOT EXISTS "${HAND_EYE_CORE_SOLVER_CPP}" OR
   NOT EXISTS "${HAND_EYE_CORE_SOLVER_HPP}")
  message(FATAL_ERROR
    "The generated ${HAND_EYE_CORE_SOLVER_NAME}.cpp/.hpp files are missing. "
    "Build the top-level generate_handeye_axxb_solver target first.")
endif()

add_executable(test_handeye_axxb
  test_handeye_axxb.cpp
  handeye_axxb.cpp
  "${HAND_EYE_CORE_SOLVER_CPP}")

target_compile_definitions(test_handeye_axxb
  PRIVATE
    POLYJAM_HAND_EYE_PAIRS=${POLYJAM_HAND_EYE_PAIRS}
    POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

target_include_directories(test_handeye_axxb
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${HAND_EYE_CORE_SOLVER_DIR}")

target_link_libraries(test_handeye_axxb PRIVATE ${POLYJAM_EIGEN_TARGET})
