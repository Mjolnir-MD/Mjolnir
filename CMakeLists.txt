cmake_minimum_required(VERSION 2.8)
enable_testing()
project(mjolnir)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extlib)

option(DEBUG "dump debug information" OFF)
if(DEBUG)
    add_definitions(-DMJOLNIR_DEBUG)
else(DEBUG)
    add_definitions(-DMJOLNIR_NO_DEBUG)
endif(DEBUG)


add_subdirectory(src)
add_subdirectory(test)
