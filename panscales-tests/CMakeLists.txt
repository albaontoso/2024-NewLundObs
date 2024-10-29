# <insert description here>

cmake_minimum_required(VERSION 3.10)

# set the project name
project(PanscalesLund VERSION 0.1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O2 -g -Wall")

#----------------------------------------------------------------------
# precision_type
option(PSDDREAL "build the code using dd_real for floating-point opertations" OFF)
option(PSQDREAL "build the code using qd_real for floating-point opertations" OFF)
option(PSQUAD   "build the code using Boost's flaot for floating-point opertations"  OFF)
option(PSDOUBLEEXP "build the code using the double_exp float type"  OFF)
option(PSMPFR4096  "build the code using the mpfr::real<4096> type"  OFF)

#----------------------------------------------------------------------
# define a parent dir so we can have educated guesses for most of the
# dependences:
#   PARENT_DIR  path to the main 2020-eeshower directory [../../]
#
# Note: relative paths are all given wrt to the sourcedir
set(PARENT_DIR ../../panscales/panscales-main/2020-eeshower CACHE PATH "Path to the main (2020-eeshower) directory")
message ("-- Looking for dependences in parent directory ${PARENT_DIR}")

# set the directory to pythia and its includes
set(PYTHIA_DIR ${PARENT_DIR}/pythia-interface/pythia CACHE PATH "Path to the pythia directory") 
message ("-- Looking for Pythia8 in directory ${PYTHIA_DIR}")

#----------------------------------------------------------------------
# find GSL (inludes and lib)
find_package(GSL REQUIRED)

#----------------------------------------------------------------------
# find CmdLine (inludes and lib)
#   CMDLINE_DIR  path to CmdLine (includes and lib)  [${PARENT_DIR}/helpers/CmdLine]

# find_path( CMDLINE_DIR
#   NAMES CmdLine.hh
#   HINTS ${PARENT_DIR}/helpers/CmdLine
#   )
# if (CMDLINE_DIR)
#   message ("-- CmdLine in ${CMDLINE_DIR}")
# else()
#   message(FATAL_ERROR "CmdLine not found")
# endif()
# 
# find_library( CMDLINE_LIBRARIES 
#   NAMES CmdLine
#   HINTS ${CMDLINE_DIR}
# )
# if (NOT CMDLINE_LIBRARIES)
#   message(FATAL_ERROR "CmdLine library not found")
# endif()

#----------------------------------------------------------------------
# find the AnalysisTools (includes)
#   ANALYSISTOOLS_DIR  path to AnalysisTools (includes and lib)  [${PARENT_DIR}/helpers/AnalysisTools]

find_path( ANALYSISTOOLS_DIR
  NAMES GSLRandom.hh
  #HINTS ${PARENT_DIR}/submodules/AnalysisTools AnalysisTools 
  HINTS ${PARENT_DIR}/submodules/AnalysisTools
)
if (ANALYSISTOOLS_DIR)
  message ("-- AnalysisTools in ${ANALYSISTOOLS_DIR}")
else()
  message(FATAL_ERROR "AnalysisTools not found")
endif()

#----------------------------------------------------------------------
# find the PanScales directories
#
# In the long run we want to support multiple precision types and
# out-of-source buils so weneed to introduce additional variables
#
#   PANSCALES_SOURCE_DIR  where most includes will be found [${PARENT_DIR}/shower-code]
#   PANSCALES_BUILD_DIR   where config.hh and the library are found [PANSCALES_SOURCE_DIR]

# find the source dir
find_path( PANSCALES_SOURCE_DIR 
  NAMES Type.hh
  HINTS ${PARENT_DIR}/shower-code
)
if (PANSCALES_SOURCE_DIR)
  message ("-- PanScales sources (includes) in ${PANSCALES_SOURCE_DIR}")
else()
  message(FATAL_ERROR "PanScales sources not found")
endif()

if (PANSCALES_SOURCE_DIR)
  message ("-- PanScales sources (includes) in ${PANSCALES_SOURCE_DIR}")
else()
  message(FATAL_ERROR "PanScales sources not found")
endif()

# find the include dirs for the config.hh file
find_path( PANSCALES_BUILD_DIR
  NAMES config.hh
  HINTS ${PANSCALES_SOURCE_DIR}/build-double
)
if (PANSCALES_BUILD_DIR)
  message ("-- PanScales build (config.hh) in ${PANSCALES_BUILD_DIR}")
else()
  message(FATAL_ERROR "PanScales build not found")
endif()

#----------------------------------------------------------------------
# find QD (inludes and lib) and MPFR
#
# We'll use the FindQD.cmake module in the PanScales shower dir
set(CMAKE_MODULE_PATH ${PANSCALES_SOURCE_DIR})
if (PSDDREAL OR PSQDREAL)
  find_package(QD REQUIRED)
endif()
if (PSMPFR4096)
  find_package(MPFR 2.3.0)  
endif()

#----------------------------------------------------------------------
# find the PanScales libraries 
#
# we want to do this after handling qd (maybe overkill?)

# find the PanScales library
find_library( PANSCALES_LIBRARIES 
  NAMES panscales
  HINTS ${PANSCALES_BUILD_DIR}
)
if (PANSCALES_LIBRARIES)
  message("-- PanScales library ${PANSCALES_LIBRARIES}")
else()
  message(FATAL_ERROR "PanScales library not found")
endif()

#-------------------------------------------------------------------------
# PYTHIA (if available)
option(WITH_PYTHIA "enables PYTHIA (if found)" OFF)
if (WITH_PYTHIA)
  #----------------------------------------------------
  # find the PanScalesPythia library (+Pythia itself)
  set(PANSCALESPYTHIA_SOURCE_DIR ${PANSCALES_SOURCE_DIR}/../pythia-interface)
  find_library( PANSCALESPYTHIA_LIBRARIES
    NAMES panscalespythia
    HINTS ${PANSCALESPYTHIA_SOURCE_DIR}/build-double ${PANSCALESPYTHIA_BUILD_DIR}
  )
  if (PANSCALESPYTHIA_LIBRARIES)
    message("-- PanScalesPythia library ${PANSCALESPYTHIA_LIBRARIES}")
  else()
    message(FATAL_ERROR "PanScalesPythia library not found")
  endif()
  set(PANSCALESPYTHIA_INCLUDES ${PANSCALESPYTHIA_SOURCE_DIR} ${PYTHIA_DIR}/include)
  message("-- PanScalesPythia includes ${PANSCALESPYTHIA_INCLUDES}")

  find_library(PYTHIA_LIBRARIES
    NAMES pythia8
    HINTS ${PYTHIA_DIR}/lib
  )
endif()

#-------------------------------------------------------------------------
# add the executable
add_executable(lund-analysis lund-analysis.cc)
target_include_directories(lund-analysis PUBLIC ${GSL_INCLUDE_DIR} ${QD_INCLUDE_DIR}  ${MPFR_INCLUDES} ${ANALYSISTOOLS_DIR} ${PANSCALES_SOURCE_DIR}  ${PANSCALES_BUILD_DIR})
target_link_libraries(lund-analysis ${PANSCALES_LIBRARIES} ${GSL_LIBRARIES} ${QD_LIBRARIES} ${MPFR_LIBRARIES})

