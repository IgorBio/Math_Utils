cmake_minimum_required(VERSION 3.13)
project(MathUtilsTests)

set(CMAKE_C_STANDARD 11)

include(ExternalProject)

set(CHECK_INSTALL_DIR ${CMAKE_BINARY_DIR}/extern)

ExternalProject_Add(check_project
    PREFIX ${CHECK_INSTALL_DIR}
    GIT_REPOSITORY https://github.com/libcheck/check.git
    GIT_TAG tags/0.15.2
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

set(CHECK_DIR ${CHECK_INSTALL_DIR}/src/check_project-build)
link_directories(${CHECK_DIR}/lib)
include_directories(${CHECK_DIR}/include)

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/../
)

add_executable(
    ${PROJECT_NAME}
    ${PROJECT_SOURCE_DIR}/tests.c
    ${PROJECT_SOURCE_DIR}/../math_utils.c
)

target_compile_options(
        ${PROJECT_NAME}
        PUBLIC
        -Wall
        -Werror
        -Wextra
        -Wpedantic
        -std=c11
        -D_GNU_SOURCE
)

target_link_libraries(${PROJECT_NAME} check)
add_dependencies(${PROJECT_NAME} check_project)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")

enable_testing()
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

add_custom_target(coverage
    COMMAND ${CMAKE_MAKE_PROGRAM} test
    COMMAND lcov --directory . --capture --output-file coverage.info --ignore-errors version
    COMMAND lcov --remove coverage.info '/usr/*' --output-file coverage.info
    COMMAND genhtml coverage.info --output-directory coverage
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating code coverage report"
)