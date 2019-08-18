# Example of command if a lib is located in a custom directory:
# cmake -DBOOST_ROOT:PATHNAME=~/libs/Boost_x11 -DBoost_LIBRARY_DIRS:FILEPATH=~/libs/Boost_x11/lib ..

CMAKE_MINIMUM_REQUIRED(VERSION 3.8.2)

ENABLE_LANGUAGE(CUDA)
# ENABLE_LANGUAGE(CXX)

PROJECT(RD_CHALLENGE)

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
SET(CMAKE_CXX_FLAGS "-O3 -Wall")

# Add the C++14 flags for this project
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED 14)

IF(POLICY CMP0042)
    CMAKE_POLICY(SET CMP0042 NEW) # MACOSX_RPATH
ENDIF()

OPTION(THREAD       "Build binaries with multi-threading enabled" ON)
OPTION(GPU          "Build binaries with GPU enabled" ON)
OPTION(DOC          "Setup the doxygen documentation (must run \"make doc\" afterwards to actually build the doc)" OFF)
OPTION(TESTS        "Setup the test suite" OFF)
OPTION(PROFILING    "Enable profiling" OFF)

# Get the last tag committed
# Another command can get the current branch: git rev-parse --abbrev-ref HEAD
# EXECUTE_PROCESS(
#     COMMAND git describe --abbrev=0 --tags --match=[0-9].[0-9].[0-9]
#     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#     OUTPUT_VARIABLE GIT_BRANCH
#     OUTPUT_STRIP_TRAILING_WHITESPACE
# )

# IF(THREAD)
#     ADD_DEFINITIONS(-DTHREAD)
# ENDIF(THREAD)

FIND_PACKAGE(Boost 1.65.0 REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(SYSTEM ${Boost_LIBRARY_DIRS})

FIND_PACKAGE(CUDA 10.1 REQUIRED)
# LIST(APPEND CUDA_NVCC_FLAGS -Xcompiler -fopenmp)

ADD_SUBDIRECTORY(src)

IF( DOC )
    ADD_SUBDIRECTORY(doc)
ENDIF( DOC )

IF( TESTS )
    ENABLE_TESTING()
    ADD_SUBDIRECTORY(tests)
ENDIF( TESTS )
