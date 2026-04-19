cmake_minimum_required(VERSION 3.16)

project(polyjam LANGUAGES CXX)

set(POLYJAM_VERSION_MAJOR 2)
set(POLYJAM_VERSION_MINOR 0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_compile_options(-w -O3)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w -O3")

option(POLYJAM_USE_OPENCV "Enable OpenCV visualization support." OFF)
option(POLYJAM_ENABLE_PARALLEL_PRUNING "Use divide-and-conquer row-basis pruning during solver generation." ON)
set(POLYJAM_TBB_ROOT "${PROJECT_SOURCE_DIR}/../third_party/oneTBB-2019_U9" CACHE PATH "Vendored oneTBB root used for Polyjam parallel pruning.")
set(POLYJAM_PARALLEL_PRUNING_THREADS 0 CACHE STRING "Maximum oneTBB parallelism for template row-basis pruning; 0 means use oneTBB default.")
set(POLYJAM_PARALLEL_PRUNING_MIN_BLOCK_ROWS 64 CACHE STRING "Minimum rows per pruning block before splitting work across oneTBB tasks.")

if(NOT DEFINED MACAULAYCOMMAND)
  set(MACAULAYCOMMAND "M2")
endif()

set(WORKSPACEPATH "${PROJECT_SOURCE_DIR}/../polyjam_workspace/")
set(SOLVERPATH "${PROJECT_SOURCE_DIR}/../polyjam_solvers/")

file(GLOB_RECURSE POLYJAM_SOURCE_FILES CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/src/*.cpp")

file(GLOB_RECURSE POLYJAM_HEADER_FILES CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/include/*.hpp")

add_library(polyjam SHARED ${POLYJAM_SOURCE_FILES} ${POLYJAM_HEADER_FILES})

if(POLYJAM_ENABLE_PARALLEL_PRUNING)
  set(POLYJAM_ENABLE_PARALLEL_PRUNING_VALUE 1)
else()
  set(POLYJAM_ENABLE_PARALLEL_PRUNING_VALUE 0)
endif()

include("${PROJECT_SOURCE_DIR}/../cmake/PolyjamVendoredTBB.cmake")

target_include_directories(polyjam PUBLIC ${PROJECT_SOURCE_DIR}/include)

target_compile_definitions(polyjam
  PUBLIC
    "MACAULAYCOMMAND=\"${MACAULAYCOMMAND}\""
    "WORKSPACEPATH=\"${WORKSPACEPATH}\""
    "SOLVERPATH=\"${SOLVERPATH}\""
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
