
# detect the environment
# CUDA_FOUND, MKL_FOUND, WIN32, UNIX, APPLE, PYTHONLIBS_FOUND

# read in what will be tested for any environment
FILE(GLOB TEST_CASE_ALWAYS RELATIVE ${CMAKE_SOURCE_DIR}/test/integration/  ./cases/*.cfg  )
message("Test cases found : ${TEST_CASE_ALWAYS}")

# according to the building environment, generate the list of available test cases

set(TEST_CASE ${TEST_CASE_ALWAYS})

message("Test cases to be performed : ${TEST_CASE}")

# write the available test cases to a file
file(WRITE ./test_cases.txt)
foreach(test ${TEST_CASE})
    file(APPEND ./test_cases.txt "${test}\n")
endforeach()
