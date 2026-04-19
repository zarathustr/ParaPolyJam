cmake_minimum_required(VERSION 3.16)

project(polyjam_handeye_axxb
  VERSION 2.4
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Requested optimization/no-warning flags for every target in this project.
add_compile_options(-w -O3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -O3")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

option(POLYJAM_USE_OPENCV "Enable OpenCV visualization support inside Polyjam." OFF)
option(POLYJAM_ENABLE_PARALLEL_PRUNING "Use divide-and-conquer row-basis pruning during Polyjam solver generation." ON)
set(POLYJAM_TBB_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/third_party/oneTBB-2019_U9" CACHE PATH "Vendored oneTBB root used for Polyjam parallel pruning.")
set(POLYJAM_PARALLEL_PRUNING_THREADS 0 CACHE STRING "Maximum oneTBB parallelism for template row-basis pruning; 0 means use oneTBB default.")
set(POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS 64 CACHE STRING "Minimum rows per pruning block before splitting work across oneTBB tasks.")
option(POLYJAM_BUILD_HAND_EYE_GENERATOR "Build the hand-eye AX=XB solver generator." ON)
option(POLYJAM_BUILD_HAND_EYE_TEST "Build the synthetic SE(3) hand-eye AX=XB solver test executable." ON)
option(POLYJAM_BUILD_WU_Q_LAMBDA_GENERATOR "Build the W*u=(Q-lambda*I)q solver generator." ON)
option(POLYJAM_BUILD_WU_Q_LAMBDA_TEST "Build the synthetic W*u=(Q-lambda*I)q solver test executable." ON)
option(POLYJAM_BUILD_WU_Q_LAMBDA0_GENERATOR "Build the lambda-free W*u=Q*q quaternion solver generator." ON)
option(POLYJAM_BUILD_WU_Q_LAMBDA0_TEST "Build the synthetic lambda-free W*u=Q*q quaternion solver test executable." ON)
option(POLYJAM_BUILD_GXY_MONOMIAL_GENERATOR "Build the MATLAB M1/M2 monomial-system solver generator." ON)
option(POLYJAM_BUILD_GXY_MONOMIAL_TEST "Build the synthetic MATLAB M1/M2 monomial-system solver test executable." ON)

set(POLYJAM_HAND_EYE_PAIRS 2 CACHE STRING "Number of A/B motion pairs compiled into the generated AX=XB solver.")
set(POLYJAM_MACAULAY_COMMAND "M2" CACHE STRING "Console command used to run Macaulay2.")
set(POLYJAM_TEST_SEEDS 100 CACHE STRING "Default number of RNG seeds used by synthetic statistical tests.")

set(POLYJAM_WORKSPACE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/")
set(POLYJAM_SOLVER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/")

set(POLYJAM_SOURCE_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/core/Coefficient.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/core/Monomial.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/core/Poly.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/core/PolyMatrix.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/core/Term.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/fields/Q.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/fields/R.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/fields/Sym.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/fields/Zp.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/generator/CMatrix.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/generator/ExportMacaulay.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/generator/ParallelPruning.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/generator/methods.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/math/GaussJordan.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/src/polyjam.cpp")

add_library(polyjam STATIC ${POLYJAM_SOURCE_FILES})

if(POLYJAM_ENABLE_PARALLEL_PRUNING)
  set(POLYJAM_ENABLE_PARALLEL_PRUNING_VALUE 1)
else()
  set(POLYJAM_ENABLE_PARALLEL_PRUNING_VALUE 0)
endif()

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/PolyjamVendoredTBB.cmake")
target_include_directories(polyjam
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_generator/include")

target_compile_definitions(polyjam
  PUBLIC
    "MACAULAYCOMMAND=\"${POLYJAM_MACAULAY_COMMAND}\""
    "WORKSPACEPATH=\"${POLYJAM_WORKSPACE_PATH}\""
    "SOLVERPATH=\"${POLYJAM_SOLVER_PATH}\""
    POLYJAM_ENABLE_PARALLEL_PRUNING=${POLYJAM_ENABLE_PARALLEL_PRUNING_VALUE}
    POLYJAM_PARALLEL_PRUNING_THREADS=${POLYJAM_PARALLEL_PRUNING_THREADS}
    POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS=${POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS})

polyjam_configure_vendored_tbb(polyjam)

if(POLYJAM_USE_OPENCV)
  find_package(OpenCV REQUIRED)
  target_compile_definitions(polyjam PUBLIC POLYJAM_USE_OPENCV=1)
  target_include_directories(polyjam PUBLIC ${OpenCV_INCLUDE_DIRS})
  target_link_libraries(polyjam PUBLIC ${OpenCV_LIBS})
else()
  target_compile_definitions(polyjam PUBLIC POLYJAM_USE_OPENCV=0)
endif()

# Generator executables -------------------------------------------------------
if(POLYJAM_BUILD_HAND_EYE_GENERATOR)
  add_executable(polyjam_generate_handeye_axxb
    "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/handeye_axxb/handeye_axxb.cpp")
  target_link_libraries(polyjam_generate_handeye_axxb PRIVATE polyjam)
  target_compile_definitions(polyjam_generate_handeye_axxb
    PRIVATE POLYJAM_HAND_EYE_PAIRS=${POLYJAM_HAND_EYE_PAIRS})
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA_GENERATOR)
  add_executable(polyjam_generate_wu_q_lambda
    "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/wu_q_lambda/wu_q_lambda.cpp")
  target_link_libraries(polyjam_generate_wu_q_lambda PRIVATE polyjam)
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA0_GENERATOR)
  add_executable(polyjam_generate_wu_q_lambda0
    "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/wu_q_lambda0/wu_q_lambda0.cpp")
  target_link_libraries(polyjam_generate_wu_q_lambda0 PRIVATE polyjam)
endif()

if(POLYJAM_BUILD_GXY_MONOMIAL_GENERATOR)
  add_executable(polyjam_generate_gxy_monomial
    "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/gxy_monomial/gxy_monomial.cpp")
  target_link_libraries(polyjam_generate_gxy_monomial PRIVATE polyjam)
endif()

# Solver paths ----------------------------------------------------------------
set(HAND_EYE_PUBLIC_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/handeye_axxb")
set(HAND_EYE_PUBLIC_SOLVER_CPP "${HAND_EYE_PUBLIC_SOLVER_DIR}/handeye_axxb.cpp")
set(HAND_EYE_PUBLIC_SOLVER_HPP "${HAND_EYE_PUBLIC_SOLVER_DIR}/handeye_axxb.hpp")

set(HAND_EYE_CORE_SOLVER_NAME "handeye_axxb_q")
set(HAND_EYE_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/${HAND_EYE_CORE_SOLVER_NAME}")
set(HAND_EYE_CORE_SOLVER_CPP "${HAND_EYE_CORE_SOLVER_DIR}/${HAND_EYE_CORE_SOLVER_NAME}.cpp")
set(HAND_EYE_CORE_SOLVER_HPP "${HAND_EYE_CORE_SOLVER_DIR}/${HAND_EYE_CORE_SOLVER_NAME}.hpp")

set(WU_Q_LAMBDA_PUBLIC_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/wu_q_lambda")
set(WU_Q_LAMBDA_PUBLIC_SOLVER_CPP "${WU_Q_LAMBDA_PUBLIC_SOLVER_DIR}/wu_q_lambda.cpp")
set(WU_Q_LAMBDA_PUBLIC_SOLVER_HPP "${WU_Q_LAMBDA_PUBLIC_SOLVER_DIR}/wu_q_lambda.hpp")

set(WU_Q_LAMBDA_CORE_SOLVER_NAME "wu_q_lambda_core")
set(WU_Q_LAMBDA_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/${WU_Q_LAMBDA_CORE_SOLVER_NAME}")
set(WU_Q_LAMBDA_CORE_SOLVER_CPP "${WU_Q_LAMBDA_CORE_SOLVER_DIR}/${WU_Q_LAMBDA_CORE_SOLVER_NAME}.cpp")
set(WU_Q_LAMBDA_CORE_SOLVER_HPP "${WU_Q_LAMBDA_CORE_SOLVER_DIR}/${WU_Q_LAMBDA_CORE_SOLVER_NAME}.hpp")

set(WU_Q_LAMBDA0_PUBLIC_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/wu_q_lambda0")
set(WU_Q_LAMBDA0_PUBLIC_SOLVER_CPP "${WU_Q_LAMBDA0_PUBLIC_SOLVER_DIR}/wu_q_lambda0.cpp")
set(WU_Q_LAMBDA0_PUBLIC_SOLVER_HPP "${WU_Q_LAMBDA0_PUBLIC_SOLVER_DIR}/wu_q_lambda0.hpp")

set(WU_Q_LAMBDA0_CORE_SOLVER_NAME "wu_q_lambda0_core")
set(WU_Q_LAMBDA0_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/${WU_Q_LAMBDA0_CORE_SOLVER_NAME}")
set(WU_Q_LAMBDA0_CORE_SOLVER_CPP "${WU_Q_LAMBDA0_CORE_SOLVER_DIR}/${WU_Q_LAMBDA0_CORE_SOLVER_NAME}.cpp")
set(WU_Q_LAMBDA0_CORE_SOLVER_HPP "${WU_Q_LAMBDA0_CORE_SOLVER_DIR}/${WU_Q_LAMBDA0_CORE_SOLVER_NAME}.hpp")

set(GXY_MONOMIAL_PUBLIC_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/gxy_monomial")
set(GXY_MONOMIAL_PUBLIC_SOLVER_CPP "${GXY_MONOMIAL_PUBLIC_SOLVER_DIR}/gxy_monomial.cpp")
set(GXY_MONOMIAL_PUBLIC_SOLVER_HPP "${GXY_MONOMIAL_PUBLIC_SOLVER_DIR}/gxy_monomial.hpp")

set(GXY_MONOMIAL_CORE_SOLVER_NAME "gxy_monomial_core")
set(GXY_MONOMIAL_CORE_SOLVER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/${GXY_MONOMIAL_CORE_SOLVER_NAME}")
set(GXY_MONOMIAL_CORE_SOLVER_CPP "${GXY_MONOMIAL_CORE_SOLVER_DIR}/${GXY_MONOMIAL_CORE_SOLVER_NAME}.cpp")
set(GXY_MONOMIAL_CORE_SOLVER_HPP "${GXY_MONOMIAL_CORE_SOLVER_DIR}/${GXY_MONOMIAL_CORE_SOLVER_NAME}.hpp")

# Solver-generation custom targets ------------------------------------------
if(POLYJAM_BUILD_HAND_EYE_GENERATOR)
  add_custom_command(
    OUTPUT "${HAND_EYE_CORE_SOLVER_CPP}" "${HAND_EYE_CORE_SOLVER_HPP}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${HAND_EYE_CORE_SOLVER_DIR}"
    COMMAND "$<TARGET_FILE:polyjam_generate_handeye_axxb>"
    DEPENDS polyjam_generate_handeye_axxb
            "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/handeye_axxb/handeye_axxb.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating quaternion-only Groebner-basis AX=XB hand-eye core solver with ${POLYJAM_HAND_EYE_PAIRS} A/B pairs"
    VERBATIM)

  add_custom_target(generate_handeye_axxb_solver
    DEPENDS "${HAND_EYE_CORE_SOLVER_CPP}" "${HAND_EYE_CORE_SOLVER_HPP}")
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA_GENERATOR)
  add_custom_command(
    OUTPUT "${WU_Q_LAMBDA_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA_CORE_SOLVER_HPP}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${WU_Q_LAMBDA_CORE_SOLVER_DIR}"
    COMMAND "$<TARGET_FILE:polyjam_generate_wu_q_lambda>"
    DEPENDS polyjam_generate_wu_q_lambda
            "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/wu_q_lambda/wu_q_lambda.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating Groebner-basis solver for W*u=(Q-lambda*I)q with u=q kron q kron q"
    VERBATIM)

  add_custom_target(generate_wu_q_lambda_solver
    DEPENDS "${WU_Q_LAMBDA_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA_CORE_SOLVER_HPP}")
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA0_GENERATOR)
  add_custom_command(
    OUTPUT "${WU_Q_LAMBDA0_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA0_CORE_SOLVER_HPP}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${WU_Q_LAMBDA0_CORE_SOLVER_DIR}"
    COMMAND "$<TARGET_FILE:polyjam_generate_wu_q_lambda0>"
    DEPENDS polyjam_generate_wu_q_lambda0
            "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/wu_q_lambda0/wu_q_lambda0.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating lambda-free Groebner-basis solver for three rows of W*u=Q*q plus q^Tq=1"
    VERBATIM)

  add_custom_target(generate_wu_q_lambda0_solver
    DEPENDS "${WU_Q_LAMBDA0_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA0_CORE_SOLVER_HPP}")
endif()

if(POLYJAM_BUILD_GXY_MONOMIAL_GENERATOR)
  add_custom_command(
    OUTPUT "${GXY_MONOMIAL_CORE_SOLVER_CPP}" "${GXY_MONOMIAL_CORE_SOLVER_HPP}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${GXY_MONOMIAL_CORE_SOLVER_DIR}"
    COMMAND "$<TARGET_FILE:polyjam_generate_gxy_monomial>"
    DEPENDS polyjam_generate_gxy_monomial
            "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_workspace/gxy_monomial/gxy_monomial.cpp"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Generating Groebner-basis solver for the six-equation MATLAB M1/M2 monomial system"
    VERBATIM)

  add_custom_target(generate_gxy_monomial_solver
    DEPENDS "${GXY_MONOMIAL_CORE_SOLVER_CPP}" "${GXY_MONOMIAL_CORE_SOLVER_HPP}")
endif()

# Eigen is only required for generated solvers and tests, not for Polyjam itself.
set(HAVE_EIGEN_FOR_POLYJAM_TEST OFF)
if(POLYJAM_BUILD_HAND_EYE_TEST OR
   POLYJAM_BUILD_WU_Q_LAMBDA_TEST OR
   POLYJAM_BUILD_WU_Q_LAMBDA0_TEST OR
   POLYJAM_BUILD_GXY_MONOMIAL_TEST)
  find_package(Eigen3 QUIET NO_MODULE)
  if(Eigen3_FOUND)
    set(HAVE_EIGEN_FOR_POLYJAM_TEST ON)
    set(POLYJAM_EIGEN_TARGET Eigen3::Eigen)
  else()
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
    find_package(Eigen QUIET)
    if(EIGEN_FOUND)
      set(HAVE_EIGEN_FOR_POLYJAM_TEST ON)
      add_library(polyjam_eigen INTERFACE)
      target_include_directories(polyjam_eigen INTERFACE "${EIGEN_INCLUDE_DIRS}" "${EIGEN_INCLUDE_DIR}")
      set(POLYJAM_EIGEN_TARGET polyjam_eigen)
    endif()
  endif()
endif()

# Test executables ------------------------------------------------------------
if(POLYJAM_BUILD_HAND_EYE_TEST)
  if(NOT HAVE_EIGEN_FOR_POLYJAM_TEST)
    message(WARNING "Eigen was not found; test_handeye_axxb will not be built.")
  elseif(NOT POLYJAM_BUILD_HAND_EYE_GENERATOR)
    message(WARNING "POLYJAM_BUILD_HAND_EYE_GENERATOR is OFF; test_handeye_axxb will not be built.")
  else()
    set_source_files_properties("${HAND_EYE_CORE_SOLVER_CPP}" "${HAND_EYE_CORE_SOLVER_HPP}"
      PROPERTIES GENERATED TRUE)

    add_executable(test_handeye_axxb
      "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/handeye_axxb/test_handeye_axxb.cpp"
      "${HAND_EYE_PUBLIC_SOLVER_CPP}"
      "${HAND_EYE_CORE_SOLVER_CPP}")

    add_dependencies(test_handeye_axxb generate_handeye_axxb_solver)

    target_include_directories(test_handeye_axxb
      PRIVATE
        "${HAND_EYE_PUBLIC_SOLVER_DIR}"
        "${HAND_EYE_CORE_SOLVER_DIR}")

    target_compile_definitions(test_handeye_axxb
      PRIVATE
        POLYJAM_HAND_EYE_PAIRS=${POLYJAM_HAND_EYE_PAIRS}
        POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

    target_link_libraries(test_handeye_axxb PRIVATE ${POLYJAM_EIGEN_TARGET})

    enable_testing()
    add_test(NAME handeye_axxb COMMAND test_handeye_axxb ${POLYJAM_TEST_SEEDS})

    add_custom_target(run_handeye_axxb_test
      COMMAND test_handeye_axxb ${POLYJAM_TEST_SEEDS}
      DEPENDS test_handeye_axxb
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA_TEST)
  if(NOT HAVE_EIGEN_FOR_POLYJAM_TEST)
    message(WARNING "Eigen was not found; test_wu_q_lambda will not be built.")
  elseif(NOT POLYJAM_BUILD_WU_Q_LAMBDA_GENERATOR)
    message(WARNING "POLYJAM_BUILD_WU_Q_LAMBDA_GENERATOR is OFF; test_wu_q_lambda will not be built.")
  else()
    set_source_files_properties("${WU_Q_LAMBDA_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA_CORE_SOLVER_HPP}"
      PROPERTIES GENERATED TRUE)

    add_executable(test_wu_q_lambda
      "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/wu_q_lambda/test_wu_q_lambda.cpp"
      "${WU_Q_LAMBDA_PUBLIC_SOLVER_CPP}"
      "${WU_Q_LAMBDA_CORE_SOLVER_CPP}")

    add_dependencies(test_wu_q_lambda generate_wu_q_lambda_solver)

    target_include_directories(test_wu_q_lambda
      PRIVATE
        "${WU_Q_LAMBDA_PUBLIC_SOLVER_DIR}"
        "${WU_Q_LAMBDA_CORE_SOLVER_DIR}")

    target_compile_definitions(test_wu_q_lambda
      PRIVATE POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

    target_link_libraries(test_wu_q_lambda PRIVATE ${POLYJAM_EIGEN_TARGET})

    enable_testing()
    add_test(NAME wu_q_lambda COMMAND test_wu_q_lambda ${POLYJAM_TEST_SEEDS})

    add_custom_target(run_wu_q_lambda_test
      COMMAND test_wu_q_lambda ${POLYJAM_TEST_SEEDS}
      DEPENDS test_wu_q_lambda
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()

if(POLYJAM_BUILD_WU_Q_LAMBDA0_TEST)
  if(NOT HAVE_EIGEN_FOR_POLYJAM_TEST)
    message(WARNING "Eigen was not found; test_wu_q_lambda0 will not be built.")
  elseif(NOT POLYJAM_BUILD_WU_Q_LAMBDA0_GENERATOR)
    message(WARNING "POLYJAM_BUILD_WU_Q_LAMBDA0_GENERATOR is OFF; test_wu_q_lambda0 will not be built.")
  else()
    set_source_files_properties("${WU_Q_LAMBDA0_CORE_SOLVER_CPP}" "${WU_Q_LAMBDA0_CORE_SOLVER_HPP}"
      PROPERTIES GENERATED TRUE)

    add_executable(test_wu_q_lambda0
      "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/wu_q_lambda0/test_wu_q_lambda0.cpp"
      "${WU_Q_LAMBDA0_PUBLIC_SOLVER_CPP}"
      "${WU_Q_LAMBDA0_CORE_SOLVER_CPP}")

    add_dependencies(test_wu_q_lambda0 generate_wu_q_lambda0_solver)

    target_include_directories(test_wu_q_lambda0
      PRIVATE
        "${WU_Q_LAMBDA0_PUBLIC_SOLVER_DIR}"
        "${WU_Q_LAMBDA0_CORE_SOLVER_DIR}")

    target_compile_definitions(test_wu_q_lambda0
      PRIVATE POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

    target_link_libraries(test_wu_q_lambda0 PRIVATE ${POLYJAM_EIGEN_TARGET})

    enable_testing()
    add_test(NAME wu_q_lambda0 COMMAND test_wu_q_lambda0 ${POLYJAM_TEST_SEEDS})

    add_custom_target(run_wu_q_lambda0_test
      COMMAND test_wu_q_lambda0 ${POLYJAM_TEST_SEEDS}
      DEPENDS test_wu_q_lambda0
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()

if(POLYJAM_BUILD_GXY_MONOMIAL_TEST)
  if(NOT HAVE_EIGEN_FOR_POLYJAM_TEST)
    message(WARNING "Eigen was not found; test_gxy_monomial will not be built.")
  elseif(NOT POLYJAM_BUILD_GXY_MONOMIAL_GENERATOR)
    message(WARNING "POLYJAM_BUILD_GXY_MONOMIAL_GENERATOR is OFF; test_gxy_monomial will not be built.")
  else()
    set_source_files_properties("${GXY_MONOMIAL_CORE_SOLVER_CPP}" "${GXY_MONOMIAL_CORE_SOLVER_HPP}"
      PROPERTIES GENERATED TRUE)

    add_executable(test_gxy_monomial
      "${CMAKE_CURRENT_SOURCE_DIR}/polyjam_solvers/gxy_monomial/test_gxy_monomial.cpp"
      "${GXY_MONOMIAL_PUBLIC_SOLVER_CPP}"
      "${GXY_MONOMIAL_CORE_SOLVER_CPP}")

    add_dependencies(test_gxy_monomial generate_gxy_monomial_solver)

    target_include_directories(test_gxy_monomial
      PRIVATE
        "${GXY_MONOMIAL_PUBLIC_SOLVER_DIR}"
        "${GXY_MONOMIAL_CORE_SOLVER_DIR}")

    target_compile_definitions(test_gxy_monomial
      PRIVATE POLYJAM_TEST_SEEDS=${POLYJAM_TEST_SEEDS})

    target_link_libraries(test_gxy_monomial PRIVATE ${POLYJAM_EIGEN_TARGET})

    enable_testing()
    add_test(NAME gxy_monomial COMMAND test_gxy_monomial ${POLYJAM_TEST_SEEDS})

    add_custom_target(run_gxy_monomial_test
      COMMAND test_gxy_monomial ${POLYJAM_TEST_SEEDS}
      DEPENDS test_gxy_monomial
      WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()
