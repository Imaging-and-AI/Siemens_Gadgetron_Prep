cmake_minimum_required(VERSION 3.10)

project(siemens_gadgetron)

# --- Provide good defaults for searching for packages (i.e. ismrmrd)
if(EXISTS $ENV{CMAKE_PREFIX_PATH})
  list(APPEND CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif()
if(EXISTS $ENV{ISMRMRD_HOME})
  list(APPEND CMAKE_PREFIX_PATH $ENV{ISMRMRD_HOME})
endif()
if (CMAKE_PREFIX_PATH)
    list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH)
endif ()
message(STATUS "Looking for packages in : ${CMAKE_PREFIX_PATH}")
# ---

# CMake 3.0 ``MACOSX_RPATH`` is enabled by default.
if (APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()
set(CMAKE_CXX_STANDARD 17)
# Require C++11 support
include(CheckCXXCompilerFlag)

# check the compiler version
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (NOT CMAKE_CXX_COMPILER_VERSION) #Compiler version is not set on Ubuntu 12.02 (gcc 4.6)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
  endif()
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    message("CXX COMPILER VERSION: ${CMAKE_CXX_COMPILER_ID} : ${CMAKE_CXX_COMPILER_VERSION}")
    message(FATAL_ERROR "Gadgetron requires GCC version >= 4.8")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
    message(FATAL_ERROR "Gadgetron requires Clang version >= 3.5")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
    message(FATAL_ERROR "Gadgetron requires MSVC 2013 or above ")
  endif()
else()
  message(WARNING "Unsupported compiler!")
endif()

find_package(Gadgetron CONFIG REQUIRED)
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

# build options for OpenMP support
find_package(OpenMP)
OPTION(USE_OPENMP "Use OpenMP" On)
if (OPENMP_FOUND)
    if(USE_OPENMP) 
        message("OpenMP multithreading enabled")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        ADD_DEFINITIONS(-DUSE_OMP)
    else (USE_OPENMP)
        message("OpenMP multithreading is supported, but disabled")
    endif(USE_OPENMP) 
else (OPENMP_FOUND)
  message("OpenMP multithreading not supported")
endif (OPENMP_FOUND)

if (WIN32)
    ADD_DEFINITIONS(-DWIN32 -D_WIN32 -D_WINDOWS)
#    ADD_DEFINITIONS(-DUNICODE -D_UNICODE)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
    ADD_DEFINITIONS(-D_VARIADIC_MAX=10) #to fix compiler limitations in Visual Studio Express
    ADD_DEFINITIONS("/wd4251") #disable warnings, 4251: needs to have dll-interface to be used by clients
    ADD_DEFINITIONS("/wd4344") #disable warnings, 4344: behavior change: use of explicit template arguments
    ADD_DEFINITIONS("/wd4996") #disable warnings, 4996: the POSIX name for this item is deprecated. Instead, use the ISO C++ conformant name
    if ( HAS_64_BIT )
        ADD_DEFINITIONS(-DWIN64 -D_WIN64)
    endif ( HAS_64_BIT )
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /MP /arch:AVX")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /bigobj")
    SET (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} /INCREMENTAL:NO")
    SET (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} /INCREMENTAL:NO")
    SET (CMAKE_EXE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    SET (CMAKE_SHARED_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    SET (CMAKE_STATIC_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    SET (CMAKE_MODULE_LINKER_FLAGS_DEBUG "/debug /INCREMENTAL:NO")
    # The two flags below is to fix Windows problems in relation to multiple defined operators new/delete and some constructors that are defined in our headers
    #SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE") 
    #SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /FORCE:MULTIPLE") 
else (WIN32)
    if (UNIX)
        if (APPLE)
            set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -stdlib=libc++")
        else ()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=return-type")
        endif ()
    endif ()
endif (WIN32)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${GADGETRON_HOME}/lib/cmake/ISMRMRD)

# whether to suppress compilation warnings
OPTION(BUILD_SUPPRESS_WARNINGS "Build package while suppressing warnings" Off)
if (BUILD_SUPPRESS_WARNINGS)
    if (WIN32)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
    elseif (WIN32)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
    endif (WIN32)
endif (BUILD_SUPPRESS_WARNINGS)

# whether to suppress compilation warnings
OPTION(BUILD_FOR_TIMING "Build package while timing purposes" Off)
if (BUILD_FOR_TIMING)
    ADD_DEFINITIONS(-DBUILD_WITH_TIMING)
endif (BUILD_FOR_TIMING)


if(WIN32)
  add_definitions( -DBOOST_ALL_NO_LIB )
  add_definitions( -DBOOST_ALL_DYN_LINK )
endif(WIN32)



IF (WIN32)
    set(CMAKE_DEBUG_SUFFIX d CACHE STRING "the debug suffix")
ENDIF(WIN32)

###############################################################
###############################################################

add_subdirectory(gadgets)
add_subdirectory(toolboxes)
add_subdirectory(test)
