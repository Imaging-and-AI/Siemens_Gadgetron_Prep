cmake_minimum_required(VERSION 3.17...3.22 FATAL_ERROR)

if (${CMAKE_VERSION} VERSION_LESS 3.17)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else ()
    cmake_policy(VERSION 3.17)
endif()

project(siemens_gadgetron LANGUAGES CXX C)


# --- Provide good defaults for searching for packages (i.e. ismrmrd)
find_program(CCACHE_FOUND ccache)
if (CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif (CCACHE_FOUND)

if (EXISTS $ENV{CMAKE_PREFIX_PATH})
    list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif ()
if (EXISTS $ENV{ISMRMRD_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ISMRMRD_HOME})
endif ()
if (EXISTS $ENV{TORCH_HOME})
    list(APPEND CMAKE_PREFIX_PATH $ENV{TORCH_HOME})
endif ()
if (CMAKE_PREFIX_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
endif ()
message(STATUS "Looking for packages in : ${CMAKE_PREFIX_PATH}")
# ---

# check the compiler version
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (NOT CMAKE_CXX_COMPILER_VERSION) #Compiler version is not set on Ubuntu 12.02 (gcc 4.6)
        execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
    endif ()
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 6.0)
        message("CXX COMPILER VERSION: ${CMAKE_CXX_COMPILER_ID} : ${CMAKE_CXX_COMPILER_VERSION}")
        message(FATAL_ERROR "Gadgetron requires GCC version >= 6.0")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        message(FATAL_ERROR "Gadgetron requires Clang version >= 3.5")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        message(FATAL_ERROR "Gadgetron requires MSVC 2013 or above ")
    endif ()
else ()
    message(WARNING "Unsupported compiler!")
endif ()

######################################################

set(CMAKE_MODULE_PATH ${GADGETRON_CMAKE_MODULES} ${ISMRMRD_CMAKE_MODULES} ${CMAKE_SOURCE_DIR}/cmake $ENV{CONDA_PREFIX}/share/gadgetron/cmake)

######################################################
#Gadgetron
find_package(Gadgetron REQUIRED)
find_package(ISMRMRD REQUIRED)

message("GADGETRON_HOME is ${GADGETRON_HOME}")
list(APPEND CMAKE_MODULE_PATH ${GADGETRON_HOME}/lib/cmake/ISMRMRD ${GADGETRON_HOME}/share/gadgetron/cmake)
message("CMAKE_MODULE_PATH is " ${CMAKE_MODULE_PATH})

######################################################

set(GADGETRON_VERSION_MAJOR 4)
set(GADGETRON_VERSION_MINOR 6)
set(GADGETRON_VERSION_PATCH 0)
set(GADGETRON_VERSION_STRING ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR}.${GADGETRON_VERSION_PATCH})
set(GADGETRON_SOVERSION ${GADGETRON_VERSION_MAJOR}.${GADGETRON_VERSION_MINOR})

# Compatibility variables from FindGadgetron.cmake
set(GADGETRON_INCLUDE_DIR         ${GADGETRON_HOME}/include/gadgetron)
set(GADGETRON_LIB_DIR             ${GADGETRON_HOME}/lib)
set(GADGETRON_INSTALL_CONFIG_PATH ${GADGETRON_HOME}/share/gadgetron/config)
set(GADGETRON_VERSION_MAJOR       ${Gadgetron_VERSION_MAJOR})
set(GADGETRON_VERSION_MINOR       ${Gadgetron_VERSION_MINOR})
set(GADGETRON_VERSION_PATCH       ${Gadgetron_VERSION_PATCH})

set(CMAKE_INSTALL_PREFIX ${GADGETRON_HOME})
message("set CMAKE_INSTALL_PREFIX to be the gadgetron home : ${CMAKE_INSTALL_PREFIX}")
LINK_DIRECTORIES(${GADGETRON_LIB_DIR})

find_package(Git)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE SOURCE_GIT_SHA1_HASH ERROR_VARIABLE GADGETRON_GIT_STDERR)
    string(STRIP "${SOURCE_GIT_SHA1_HASH}" SOURCE_GIT_SHA1_HASH)
    string(LENGTH "${SOURCE_GIT_SHA1_HASH}" SOURCE_GIT_SHA1_HASH_LEN)
    if (${SOURCE_GIT_SHA1_HASH_LEN} LESS 40)
        message(WARNING "Could not determine SHA-1 hash: ${GADGETRON_GIT_STDERR}")
        set(SOURCE_GIT_SHA1_HASH "NA")
    endif ()
else ()
    set(SOURCE_GIT_SHA1_HASH "NA")
endif ()

set(GADGETRON_INSTALL_CMAKE_PATH share/gadgetron/cmake)
set(GADGETRON_INSTALL_CONFIG_PATH share/gadgetron/config)
set(GADGETRON_INSTALL_LOGO_PATH share/gadgetron/logo)
set(GADGETRON_INSTALL_MATLAB_PATH share/gadgetron/matlab)
set(GADGETRON_INSTALL_PYTHON_MODULE_PATH share/gadgetron/python)
set(GADGETRON_INSTALL_SCHEMA_PATH schema)
set(GADGETRON_INSTALL_INCLUDE_PATH include/gadgetron)
set(GADGETRON_INSTALL_CHROOT_SCRIPTS_PATH share/gadgetron/chroot)
if (WIN32)
    set(GADGETRON_VAR_DIR $ENV{APPDATA}/gadgetron CACHE PATH "Path where Gadgetron will store calibration data, etc" FORCE)
else ()
    set(GADGETRON_VAR_DIR /var/lib/gadgetron/ CACHE PATH "Path where Gadgetron will store calibration data, etc" FORCE)
endif ()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Set the build type to Release if not specified
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# build options for 64 bits system
if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("64bit system is found")
  set( HAS_64_BIT On CACHE BOOL "64bit build")
else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
  message("32bit system is found")
  set( HAS_64_BIT Off CACHE BOOL "64bit build")
endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )

set(CMAKE_INSTALL_RPATH ".:${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# build options for OpenMP support
find_package(OpenMP)
option(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if (USE_OPENMP)
        message("OpenMP multithreading enabled")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSE_OMP)
    else ()
        message("OpenMP multithreading is supported, but disabled")
    endif ()
else ()
    message("OpenMP multithreading not supported")
endif ()

set(CMAKE_CXX_STANDARD 17)
if (WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    add_definitions(-DWIN32 -D_WIN32 -D_WINDOWS -DWIN -D_AMD64_)
    add_definitions(-Zc:__cplusplus)
    add_definitions(-DNOMINMAX)
    #add_definitions(-DUNICODE -D_UNICODE)
    add_definitions("/wd4251") #disable warnings, 4251: needs to have dll-interface to be used by clients/
    if (HAS_64_BIT)
        add_definitions(-DWIN64 -D_WIN64)
    endif ()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /bigobj /openmp:experimental /permissive- /FS /DNOGDI /we4667")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /EHsc /MP /bigobj /w1")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} /INCREMENTAL:NO")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_STATIC_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    set(CMAKE_MODULE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    # The two flags below is to fix Windows problems in relation to multiple defined operators new/delete and some constructors that are defined in our headers
    #  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE")
    #  set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE")
else ()
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    include(CheckCXXCompilerFlag)
    if (APPLE)
        if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif ()
    endif ()

    check_cxx_compiler_flag(-Werror=return-type COMPILER_SUPPORTS_WRETURNTYPE)
    if (COMPILER_SUPPORTS_WRETURNTYPE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
    endif ()
    check_cxx_compiler_flag(-Werror=terminate COMPILER_SUPPORTS_WTERMINATE)
    if (COMPILER_SUPPORTS_WTERMINATE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=terminate")
    endif ()
    check_cxx_compiler_flag(-Werror=non-virtual-dtor COMPILER_SUPPORT_WNON_VIRTUAL)
    if (COMPILER_SUPPORT_WNON_VIRTUAL)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=non-virtual-dtor")
    endif ()
    
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif ()

# whether to suppress compilation warnings
option(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    if (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    elseif (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif ()
endif ()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# uncomment these if all compiled targets are to be stored in the same directory
if (WIN32)
    SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/prod/bin CACHE STRING "Where to put the executables")
    SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/prod/bin CACHE STRING "Where to put the libraries")
endif (WIN32)

# whether to suppress compilation warnings
OPTION(BUILD_FOR_TIMING "Build package while timing purposes" Off)
if (BUILD_FOR_TIMING)
    ADD_DEFINITIONS(-DBUILD_WITH_TIMING)
endif (BUILD_FOR_TIMING)

###############################################################
#Bootstrap search for libraries (We need to find cmake modules in Gadgetron)
###############################################################

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# necessary for Windows and RHEL <=6 systems

if(WIN32)
    add_definitions( -DHAVE_SNPRINTF) #Python decides to overwrite snprintf if this is not defined, because ofcourse it does.
endif()
# We actually only use system and thread explicitly, but they require linking in date_time and chrono
if (WIN32)
    find_package(Boost 1.71.0 COMPONENTS coroutine system date_time chrono program_options filesystem timer exception REQUIRED )
else()
    find_package(Boost 1.71.0 COMPONENTS coroutine system date_time program_options filesystem timer REQUIRED )
endif()
add_definitions("-DBOOST_BIND_GLOBAL_PLACEHOLDERS") # To suppress messages in later versions of boost, in version 1.76.0 there are internal conflicts in boost
add_definitions("-DBOOST_ALLOW_DEPRECATED_HEADERS") # To suppress messages in later versions of boost, in version 1.76.0 there are internal conflicts in boost
if (NOT DEFINED Boost_VERSION_STRING)
    set(Boost_VERSION_STRING "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}")
endif()
if (Boost_FOUND)
    message("Found Boost : ${Boost_LIBRARIES}")
endif()

###############################################################

find_package(GTest)
#Add support for the default ubuntu package of gtest (which is not compiled
if (NOT GTEST_FOUND)
  find_path(GTEST_SRC_DIR src/gtest.cc HINTS /usr/src/gtest)
  find_path(GTEST_INCLUDE_DIRS gtest.h HINTS /usr/include/gtest)
  if (GTEST_SRC_DIR AND GTEST_INCLUDE_DIRS)
    MESSAGE("GTest src package found. Compiling as part of Gadgetron.")
    add_subdirectory(${GTEST_SRC_DIR} ${CMAKE_BINARY_DIR}/gtest )
    include_directories(${GTEST_INCLUDE_DIRS})
    set(GTEST_FOUND 1)
    set(GTEST_LIBRARIES gtest gtest_main)
  endif (GTEST_SRC_DIR AND GTEST_INCLUDE_DIRS)
endif (NOT GTEST_FOUND)

###############################################################
find_package(Armadillo REQUIRED )
add_library(armadillo INTERFACE)
set_property(TARGET armadillo PROPERTY INTERFACE_include_directories ${ARMADILLO_INCLUDE_DIR})
add_definitions(-DARMA_DONT_USE_WRAPPER -DARMA_USE_CXX11 -DARMA_64BIT_WORD)
install(TARGETS armadillo EXPORT gtprep-export)

###############################################################

find_package(FFTW3 REQUIRED COMPONENTS single double)
add_library(FFTW INTERFACE)
target_link_libraries(FFTW INTERFACE ${FFTW3_LIBRARIES})
set_property(TARGET FFTW PROPERTY INTERFACE_include_directories ${FFTW3_INCLUDE_DIR})
install(TARGETS FFTW EXPORT gadgetron-export)

###############################################################
find_package(ISMRMRD REQUIRED)
if (ISMRMRD_FOUND)
    message("ISMRMRD Found")
else()
    message("ISMRMRD not found")
endif()

###############################################################
find_package(Eigen3 REQUIRED)
if (Eigen3_FOUND)
    message("Eigen3 Found")
    include_directories(${EIGEN3_INCLUDE_DIR})
else()
    message("Eigen3 not found")
endif()

###############################################################
option(USE_MKL "Use MKL for numerical computing" On)
if (USE_MKL)
    if (LINUX)
        list(APPEND CMAKE_PREFIX_PATH "/opt/intel/oneapi/mkl/latest")
    endif()
    find_package(MKL REQUIRED)
    link_directories(GTBLAS BEFORE INTERFACE "${MKL_ROOT}/lib/intel64/")
    add_library(GTBLAS INTERFACE)
    set_property(TARGET GTBLAS PROPERTY INTERFACE_LINK_LIBRARIES ${MKL_LIBRARIES})
    message("MKL LIB DIR ${MKL_ROOT}/lib/intel64")

    add_definitions(-DUSE_MKL -DARMA_BLAS_LONG_LONG -DARMA_USE_BLAS -DARMA_USE_LAPACK)
else ()
    find_package(CBLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(LAPACKE REQUIRED)
    add_library(GTBLAS INTERFACE)
    set_property(TARGET GTBLAS PROPERTY INTERFACE_include_directories ${CBLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR} ${LAPACKE_INCLUDE_DIR})
    set_property(TARGET GTBLAS PROPERTY INTERFACE_LINK_LIBRARIES ${CBLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${LAPACKE_LIBRARIES})

    find_package(PkgConfig)
    pkg_check_modules(PC_lapacke lapacke)

    string(COMPARE EQUAL ${BLA_VENDOR} "OpenBLAS" USING_OPENBLAS)
    if (USING_OPENBLAS)
        message("OpenBLAS detected; setting max number of OpenBLAS pthreads to one.")
        add_definitions(-DFORCE_LIMIT_OPENBLAS_NUM_THREADS)
    endif()
endif ()

install(TARGETS GTBLAS EXPORT gtprep-export)

###############################################################

find_package(GMatlab)
if (MATLAB_FOUND)
  add_definitions(-DUSE_MATLAB)
endif ()

###############################################################

find_package(PLplot)
if (PLPLOT_FOUND)
    include_directories(${PLPLOT_INCLUDE_DIR})
    link_directories(${PLPLOT_LIB_DIR})
endif ()

###############################################################

find_package(DCMTK)
if(DCMTK_FOUND)
    message("Compiling DICOM tests")
    ADD_DEFINITIONS(-DUSE_DCMTK)
endif ()

###############################################################

find_package(HDF5 1.8 COMPONENTS C REQUIRED)
if (HDF5_FOUND)
   message("HDF5 Found")
else()
   message("HDF5 not found")
endif()

###############################################################

find_package(Python3 COMPONENTS Development NumPy REQUIRED)
find_package(Boost COMPONENTS python3${Python3_VERSION_MINOR} REQUIRED)

message("Found numpy : ${Python3_NumPy_INCLUDE_DIRS} ${Python3_NumPy_LIBRARIES} ${Python3_INCLUDE_DIRS} ${Python3_LIBRARIES}")

if (Python3_VERSION_MINOR LESS 7)
    message(FATAL_ERROR "Gadgetron Python support requires Python 3.7 or newer. Found: ${Python3_VERSION}")
endif()

message("Boost Version: ${Boost_VERSION_STRING} (with Python ${Python3_VERSION})")
message("Secret: ${Boost_PYTHON3_TARGET}")

add_library(python INTERFACE)
set_property(TARGET python PROPERTY INTERFACE_include_directories ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})
target_link_libraries(python INTERFACE Python3::Python Python3::NumPy Boost::python3${Python3_VERSION_MINOR})

install(TARGETS python EXPORT gtprep-export)

if (Boost_FOUND)
    add_library(Boost INTERFACE)
    set_property(TARGET Boost PROPERTY INTERFACE_include_directories ${BOOST_INCLUDE_DIRS})
    target_link_libraries(Boost INTERFACE ${Boost_LIBRARIES})
endif()

###############################################################

if(EXISTS $ENV{TORCH_HOME})
    find_package(Torch)
    if (TORCH_FOUND)
        add_definitions(-DUSE_TORCH)
        message("Find Pytorch c++ library at ${TORCH_LIBRARIES}")
    else()
        message("Do not find Pytorch c++ library")
    endif ()
else()
    message("Do not find TORCH_HOME")
endif()

###############################################################

IF (WIN32)
    set(CMAKE_DEBUG_SUFFIX d CACHE STRING "the debug suffix")
ENDIF(WIN32)

###############################################################
###############################################################

add_subdirectory(gadgets)
add_subdirectory(toolboxes)
add_subdirectory(test)
