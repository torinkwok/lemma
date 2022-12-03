cmake_minimum_required(VERSION 3.10)

# SetAvg project name, version and languages here.
# Version numbers are available by including "exampleConfig.h" in
# the source. See exampleConfig.h.in for some more details.
project(autodidact)

set(HOME_DIR ${CMAKE_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set default cmake build type to release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

message("CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

# The outcome of indiscriminating build types is endless confusion. I spend a couple hours on fixing a bug caused
# by the different build types alternatively put their artifacts into the same build directory.
if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/build/debug")
elseif (CMAKE_BUILD_TYPE MATCHES Release)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/build/release")
elseif (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/build/release-dbinfo")
elseif (CMAKE_BUILD_TYPE MATCHES MinSizeRel)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/build/release-min")
endif ()

# Set executable installation directory, by default to be /usr/local/bin
set(INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}")

# Set file paths
add_definitions(-DAUTODIDACT_DIR_LOG="${CMAKE_CURRENT_SOURCE_DIR}/logs")
add_definitions(-DAUTODIDACT_DIR_DATA_ABS="${CMAKE_CURRENT_SOURCE_DIR}/data/abstraction")
add_definitions(-DAUTODIDACT_DIR_DATA_STG="${CMAKE_CURRENT_SOURCE_DIR}/data/strategy")
add_definitions(-DAUTODIDACT_DIR_DATA_LAB="${CMAKE_CURRENT_SOURCE_DIR}/data/lab")
add_definitions(-DAUTODIDACT_DIR_CFG_GAME="${CMAKE_CURRENT_SOURCE_DIR}/config/games")
add_definitions(-DAUTODIDACT_DIR_CFG_ENG="${CMAKE_CURRENT_SOURCE_DIR}/config/engine")
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEV=2)
    #    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
    #    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=[sanitizer_name] [additional_options] [-g] [-OX]")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer -O1")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
    #    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
elseif (CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    add_definitions(-DDEV=1)
endif ()

if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(WARNING "please upgrade cmake")
endif ()

# Fix RPATH https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${HOME_DIR}/3rdparty/Catch2/contrib/")

if (UNIX AND NOT APPLE)
    # NOTE(kwok): To fix linker craps on Linux like this (https://stackoverflow.com/a/46811527):
    #
    #       * /lib64/libm.so.6: version `GLIBC_2.29' not found (required by release-min/agent)
    #       * /lib64/libstdc++.so.6: version `GLIBCXX_3.4.26' not found
    #
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
    set(CMAKE_SHARED_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
endif ()

# Add submodule dependencies
add_subdirectory(3rdparty/Catch2)
# Submodules that are header only
include_directories(3rdparty/cxxopts/include)
include_directories(3rdparty/spdlog/include)
include_directories(3rdparty/cereal/include)
option(THREAD_SAFE ON)

find_package(OpenSSL REQUIRED)
message(STATUS "OpenSSL Libraries Found: ${OPENSSL_FOUND}")
message(STATUS "OpenSSL Include:         ${OPENSSL_INCLUDE_DIR}")
message(STATUS "OpenSSL Libraries:       ${OPENSSL_LIBRARIES}")

# These packages are required to be installed locally, i.e. using brew, apt-get, or vspkg.
if (UNIX AND NOT APPLE)
    set(cpprestsdk_DIR /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/cmake/)
    #find_package(UUID REQUIRED uuid)
elseif (APPLE)
    #this is not safe, should add a FindTBB.cmake to automatially find and include the right directories
    #do something like this https://github.com/justusc/FindTBB
    # include_directories(/usr/local/Cellar/parallelstl/20200330/stdlib)
    # include_directories(/usr/local/Cellar/parallelstl/20200330/include)
    include_directories(${OPENSSL_INCLUDE_DIR})
    link_libraries(${OPENSSL_LIBRARIES})
endif ()
find_package(cpprestsdk REQUIRED)
find_package(Threads REQUIRED)

# Set up boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED system thread date_time regex serialization)
message(STATUS "Boost Include:      ${Boost_INCLUDE_DIR}")
message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost Libraries:    ${Boost_LIBRARIES}")

message(STATUS "C++ Flags:${CMAKE_CXX_FLAGS}")

add_subdirectory(modules/core)
add_subdirectory(modules/engine)
add_subdirectory(modules/agent)

add_subdirectory(cuckoo_lab)

# Enable "make test" command
enable_testing()
