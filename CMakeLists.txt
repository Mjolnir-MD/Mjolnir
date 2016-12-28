cmake_minimum_required(VERSION 2.8)
enable_testing()
project(mjolnir)

include_directories(${PROJECT_SOURCE_DIR})

add_subdirectory(src)
add_subdirectory(test)
