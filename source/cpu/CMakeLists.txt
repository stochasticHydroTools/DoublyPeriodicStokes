cmake_minimum_required(VERSION 3.0)
project(FCM CXX)
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)

# lib and bin will be in top source
set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")
# set runtime path for shared libraries
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.5.0")
    message(FATAL_ERROR "Insufficient gcc version. Must be 7.5.0 or greater")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(BLA_VENDOR Intel10_64_dyn)
else()
  message(FATAL_ERROR "Compiler ID must be GNU or Intel")
endif()

find_package(OpenMP REQUIRED)
find_package(LAPACK REQUIRED)

# add support for FFTW to find_package for GNU
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  configure_file(cmake.in/download_FindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
      else()
      message("CMake step for findFFTW completed (${result}).")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
          RESULT_VARIABLE result
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
  if(result)
      message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
  endif()
  
  set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)
  
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
  # now find fftw
  find_package(FFTW REQUIRED)
else()
  set(fftw3_omp "${CMAKE_INSTALL_PREFIX}/fftw_install/lib/libfftw3_omp.so")
  set(fftw3 "${CMAKE_INSTALL_PREFIX}/fftw_install/lib/libfftw3.so")
endif()

# general compiler flags
option(debug "Build with debugger" OFF)
option(fftw_wisdom "Build with FFTW wisdom enabled" ON)
option(fftw_wisdom_patient "Use fftw_patient planner" OFF)
option(use_stack "Use stack allocation for spread/interp" ON)

# setting include dir and general compile flags
# include dirs
include_directories(include)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(NOT debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Iinclude -O3 -march=native")
  elseif(debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Iinclude -O3 -g -DDEBUG -march=native")
  endif()
  include_directories(/usr/include)
  include_directories(/usr/include/openblas)
else()
  if(NOT debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Iinclude -O3 -ipo -DINTEL")
  elseif(debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Iinclude -g -DDEBUG -DINTEL")
  endif()
  include_directories(${MKLROOT}/include)
  include_directories(${CMAKE_INSTALL_PREFIX}/fftw_install/include)
endif()

if (fftw_wisdom)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_WISDOM")
  if (fftw_wisdom_patient)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW_PATIENT")
    message("Setting FFTW planner flag to FFTW_PATIENT")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_FFTW_MEASURE")
    message("Setting FFTW planner flag to FFTW_MEASURE")
  endif()
endif()

if (use_stack)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_STACK")
endif()



# setting lib src files
set(chebSRC src/Quadrature.cpp)
set(gridSRC src/Grid.cpp wrapper/GridWrapper.cpp)
set(particlesSRC src/ParticleList.cpp wrapper/ParticleListWrapper.cpp)
set(spreadInterpSRC src/SpreadInterp.cpp wrapper/SpreadInterpWrapper.cpp)
set(transformSRC src/Transform.cpp wrapper/TransformWrapper.cpp)
set(linSolveSRC src/LinearSolvers.cpp)
set(dpToolsSRC src/DPTools.cpp)
set(bcSRC wrapper/BCWrapper.cpp)


# set build output paths
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# building lib
add_library(cheb SHARED ${chebSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${chebSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
  target_link_libraries(cheb gomp) 
else()
  set_source_files_properties(${chebSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
  target_link_libraries(cheb ${OpenMP_iomp5_LIBRARY} ${OpenMP_pthread_LIBRARY})
endif()

add_library(grid SHARED ${gridSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${gridSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
  target_link_libraries(grid cheb fftw3)
else()
  set_source_files_properties(${gridSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic")
  target_link_libraries(grid cheb ${fftw3})
endif()


add_library(particles SHARED ${particlesSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${particlesSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${particlesSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(particles grid cheb)

add_library(spreadInterp SHARED ${spreadInterpSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${spreadInterpSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
else()
  set_source_files_properties(${spreadInterpSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
endif()
target_link_libraries(spreadInterp particles)

add_library(spreadInterpDerivX SHARED ${spreadInterpSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(spreadInterpDerivX PROPERTIES COMPILE_FLAGS "-DTORQUE_X -fPIC -fopenmp")
else()
  set_target_properties(spreadInterpDerivX PROPERTIES COMPILE_FLAGS "-xHost -DTORQUE_X -qopenmp -fpic")
endif()
target_link_libraries(spreadInterpDerivX particles)

add_library(spreadInterpDerivY SHARED ${spreadInterpSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(spreadInterpDerivY PROPERTIES COMPILE_FLAGS "-DTORQUE_Y -fPIC -fopenmp")
else()
  set_target_properties(spreadInterpDerivY PROPERTIES COMPILE_FLAGS "-xHost -DTORQUE_Y -qopenmp -fpic")
endif()
target_link_libraries(spreadInterpDerivY particles)

add_library(spreadInterpDerivZ SHARED ${spreadInterpSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(spreadInterpDerivZ PROPERTIES COMPILE_FLAGS "-DTORQUE_Z")
else()
  set_target_properties(spreadInterpDerivZ PROPERTIES COMPILE_FLAGS "-xHost -DTORQUE_Z -qopenmp -fpic")
endif()
target_link_libraries(spreadInterpDerivZ particles)

add_library(transform SHARED ${transformSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${transformSRC} PROPERTIES COMPILE_FLAGS "-lfftw3_omp -lfftw3 -lm -fopenmp -fPIC ")
  target_link_libraries(transform fftw3_omp fftw3 gomp)
else()
  set_source_files_properties(${transformSRC} PROPERTIES COMPILE_FLAGS "-xHost -lfftw3_omp -lfftw3 -qopenmp -fpic ")
  target_link_libraries(transform ${fftw3_omp} ${fftw3} ${OpenMP_iomp5_LIBRARY} ${OpenMP_pthread_LIBRARY})
endif()

add_library(linSolve SHARED ${linSolveSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${linSolveSRC} PROPERTIES COMPILE_FLAGS "-Winline -lm -lblas -llapacke -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -lfftw3 -fopenmp -fPIC")
  target_link_libraries(linSolve gomp openblas lapack lapacke fftw3)
else()
  set_source_files_properties(${linSolveSRC} PROPERTIES COMPILE_FLAGS "-lmkl_rt -Wl,--no-as-needed -lpthread -ldl -qopenmp -lfftw3 -fpic")
  target_link_libraries(linSolve ${BLAS_mkl_rt_LIBRARY} ${OpenMP_iomp5_LIBRARY} ${OpenMP_pthread_LIBRARY} ${fftw3})
endif()

add_library(dpTools SHARED ${dpToolsSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${dpToolsSRC} PROPERTIES COMPILE_FLAGS "-Winline -lm -lblas -llapacke -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -fopenmp -lfftw3_omp -lfftw3 -fPIC")
  target_link_libraries(dpTools linSolve fftw3_omp)
else()
  set_source_files_properties(${dpToolsSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -lfftw3 -fpic")
  target_link_libraries(dpTools linSolve)
endif()


add_library(BC SHARED ${bcSRC})
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_source_files_properties(${bcSRC} PROPERTIES COMPILE_FLAGS "-fPIC -fopenmp")
  target_link_libraries(BC gomp)
else()
  set_source_files_properties(${bcSRC} PROPERTIES COMPILE_FLAGS "-xHost -qopenmp -fpic ")
  target_link_libraries(BC ${OpenMP_iomp5_LIBRARY} ${OpenMP_pthread_LIBRARY})
endif()
# install libs
install(TARGETS spreadInterp ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS spreadInterpDerivX ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS spreadInterpDerivY ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS spreadInterpDerivZ ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS cheb ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS grid ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS particles ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS transform ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS linSolve ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS dpTools ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(TARGETS BC ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
# add env config script to top src dir
configure_file(cmake.in/config.sh.in ${CMAKE_INSTALL_PREFIX}/examples/config.sh @ONLY)
