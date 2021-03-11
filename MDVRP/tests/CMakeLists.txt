add_executable(do_tests do_tests.cpp)
target_link_libraries(do_tests PUBLIC ga gtest)
add_test(NAME TESTS COMMAND do_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})