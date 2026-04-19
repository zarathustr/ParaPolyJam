# Helper for using the bundled Intel oneTBB 2019 U9 source package.
# It builds libtbb from the vendored source tree and links Polyjam against it,
# so parallel pruning does not require OpenMP or a system TBB installation.

include_guard(GLOBAL)
include(ExternalProject)
include(ProcessorCount)

option(POLYJAM_USE_BUNDLED_TBB "Build and link the bundled oneTBB source for parallel pruning." ON)
option(POLYJAM_REQUIRE_TBB "Fail configuration when parallel pruning is enabled but TBB cannot be configured." ON)
set(POLYJAM_TBB_ROOT "" CACHE PATH "Path to an Intel TBB/oneTBB source tree. Empty uses third_party/oneTBB-2019_U9.")
set(POLYJAM_TBB_BUILD_JOBS 2 CACHE STRING "Number of jobs used to build bundled oneTBB; set higher for faster local builds.")

function(polyjam_configure_tbb out_target)
  if(TARGET polyjam_tbb)
    set(${out_target} polyjam_tbb PARENT_SCOPE)
    return()
  endif()

  if(NOT POLYJAM_USE_BUNDLED_TBB)
    find_package(TBB QUIET COMPONENTS tbb)
    if(TBB_FOUND)
      if(TARGET TBB::tbb)
        set(${out_target} TBB::tbb PARENT_SCOPE)
      elseif(DEFINED TBB_IMPORTED_TARGETS AND TBB_IMPORTED_TARGETS)
        list(GET TBB_IMPORTED_TARGETS 0 _first_tbb_target)
        set(${out_target} ${_first_tbb_target} PARENT_SCOPE)
      else()
        add_library(polyjam_tbb INTERFACE IMPORTED GLOBAL)
        target_include_directories(polyjam_tbb INTERFACE ${TBB_INCLUDE_DIRS})
        target_link_libraries(polyjam_tbb INTERFACE ${TBB_LIBRARIES})
        set(${out_target} polyjam_tbb PARENT_SCOPE)
      endif()
      return()
    endif()

    if(POLYJAM_REQUIRE_TBB)
      message(FATAL_ERROR "POLYJAM_USE_BUNDLED_TBB is OFF, but a system TBB package was not found.")
    endif()
    set(${out_target} "" PARENT_SCOPE)
    return()
  endif()

  get_filename_component(_polyjam_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
  if(POLYJAM_TBB_ROOT STREQUAL "")
    set(_tbb_root "${_polyjam_root}/third_party/oneTBB-2019_U9")
  else()
    set(_tbb_root "${POLYJAM_TBB_ROOT}")
  endif()
  get_filename_component(_tbb_root "${_tbb_root}" ABSOLUTE)

  if(NOT EXISTS "${_tbb_root}/Makefile" OR NOT EXISTS "${_tbb_root}/include/tbb/parallel_for.h")
    if(POLYJAM_REQUIRE_TBB)
      message(FATAL_ERROR "Bundled oneTBB source tree was not found or is incomplete: ${_tbb_root}")
    endif()
    set(${out_target} "" PARENT_SCOPE)
    return()
  endif()

  if(WIN32)
    find_program(POLYJAM_TBB_MAKE_TOOL NAMES gmake make)
  else()
    find_program(POLYJAM_TBB_MAKE_TOOL NAMES make gmake)
  endif()
  if(NOT POLYJAM_TBB_MAKE_TOOL)
    if(POLYJAM_REQUIRE_TBB)
      message(FATAL_ERROR "A make-compatible tool is required to build the bundled oneTBB source.")
    endif()
    set(${out_target} "" PARENT_SCOPE)
    return()
  endif()

  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(_tbb_compiler gcc)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(_tbb_compiler clang)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    if(WIN32)
      set(_tbb_compiler icl)
    else()
      set(_tbb_compiler icc)
    endif()
  elseif(MSVC)
    set(_tbb_compiler cl)
  else()
    set(_tbb_compiler gcc)
  endif()

  set(_jobs ${POLYJAM_TBB_BUILD_JOBS})
  if(_jobs LESS_EQUAL 0)
    ProcessorCount(_detected_jobs)
    if(_detected_jobs EQUAL 0)
      set(_detected_jobs 1)
    endif()
    set(_jobs ${_detected_jobs})
  endif()

  set(_tbb_build_dir "${CMAKE_BINARY_DIR}/polyjam_onetbb_build")
  set(_tbb_build_prefix "polyjam_onetbb")
  set(_tbb_release_dir "${_tbb_build_dir}/${_tbb_build_prefix}_release")

  if(WIN32)
    set(_tbb_location "${_tbb_release_dir}/tbb.dll")
    set(_tbb_import_library "${_tbb_release_dir}/tbb.lib")
  elseif(APPLE)
    set(_tbb_location "${_tbb_release_dir}/libtbb.dylib")
  else()
    set(_tbb_location "${_tbb_release_dir}/libtbb.so.2")
  endif()

  ExternalProject_Add(polyjam_onetbb_external
    SOURCE_DIR "${_tbb_root}"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "${POLYJAM_TBB_MAKE_TOOL}" -C "${_tbb_root}" -j${_jobs}
      compiler=${_tbb_compiler}
      stdver=c++11
      "CXXFLAGS=-fpermissive -w"
      tbb_build_dir=${_tbb_build_dir}
      tbb_build_prefix=${_tbb_build_prefix}
      tbb
    BUILD_BYPRODUCTS "${_tbb_location}"
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE OFF)

  add_library(polyjam_tbb SHARED IMPORTED GLOBAL)
  add_dependencies(polyjam_tbb polyjam_onetbb_external)
  set_target_properties(polyjam_tbb PROPERTIES
    IMPORTED_LOCATION "${_tbb_location}"
    INTERFACE_INCLUDE_DIRECTORIES "${_tbb_root}/include")

  if(WIN32 AND DEFINED _tbb_import_library)
    set_target_properties(polyjam_tbb PROPERTIES IMPORTED_IMPLIB "${_tbb_import_library}")
  endif()

  if(UNIX AND NOT APPLE)
    set_target_properties(polyjam_tbb PROPERTIES INTERFACE_LINK_LIBRARIES "pthread;dl")
  elseif(APPLE)
    set_target_properties(polyjam_tbb PROPERTIES INTERFACE_LINK_LIBRARIES "pthread")
  endif()

  message(STATUS "Polyjam parallel row-basis pruning will use bundled oneTBB: ${_tbb_root}")
  set(${out_target} polyjam_tbb PARENT_SCOPE)
endfunction()

function(polyjam_configure_vendored_tbb target_name)
  if(NOT TARGET ${target_name})
    message(FATAL_ERROR "polyjam_configure_vendored_tbb called for unknown target '${target_name}'")
  endif()

  if(NOT POLYJAM_ENABLE_PARALLEL_PRUNING)
    target_compile_definitions(${target_name} PUBLIC POLYJAM_HAVE_TBB=0)
    return()
  endif()

  polyjam_configure_tbb(_polyjam_tbb_target)
  if(_polyjam_tbb_target)
    target_link_libraries(${target_name} PUBLIC ${_polyjam_tbb_target})
    target_compile_definitions(${target_name} PUBLIC POLYJAM_HAVE_TBB=1)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${target_name} PUBLIC -fpermissive)
    endif()
    add_dependencies(${target_name} ${_polyjam_tbb_target})
    message(STATUS "Polyjam parallel row-basis pruning enabled with oneTBB.")
  elseif(POLYJAM_REQUIRE_TBB)
    message(FATAL_ERROR "POLYJAM_ENABLE_PARALLEL_PRUNING is ON but oneTBB could not be configured.")
  else()
    target_compile_definitions(${target_name} PUBLIC POLYJAM_HAVE_TBB=0)
    message(WARNING "oneTBB was not configured; Polyjam row-basis pruning will use the serial fallback.")
  endif()
endfunction()
