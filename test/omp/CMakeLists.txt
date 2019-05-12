set(TEST_NAMES
    test_omp_sort
    test_omp_random_number_generator
    )

add_definitions("-O2")

if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    message(STATUS "boost_unit_test_framework precompiled library found -> ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
    add_definitions(-DBOOST_TEST_DYN_LINK)
    add_definitions(-DUNITTEST_FRAMEWORK_LIBRARY_EXIST)
endif()

foreach(TEST_NAME ${TEST_NAMES})
    add_executable(${TEST_NAME} ${TEST_NAME}.cpp)
    target_link_libraries(${TEST_NAME} ${OpenMP_CXX_LIBRARIES})
    if (Boost_UNIT_TEST_FRAMEWORK_FOUND)
        target_link_libraries(${TEST_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    endif (Boost_UNIT_TEST_FRAMEWORK_FOUND)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME} --log_level=warning WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endforeach(TEST_NAME)
