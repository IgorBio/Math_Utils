cmake_minimum_required(VERSION 3.13)
project(MathUtils VERSION 1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.h
)

set(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.c
)

set(BUILD_FILES
  ${HEADERS}
  ${SOURCES}
)

add_library(${PROJECT_NAME}
  STATIC
  ${BUILD_FILES}
)

set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${CMAKE_PROJECT_NAME})
set_target_properties(${PROJECT_NAME} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

target_compile_options(
  ${PROJECT_NAME}
  PRIVATE
  -Wall
  -Werror
  -Wextra
  -Wpedantic
  -fsanitize=address
)

target_link_libraries(
  ${PROJECT_NAME}
  PRIVATE
  -fsanitize=address
)

find_program(CPPCHECK cppcheck)
if(CPPCHECK)
  message(STATUS "cppcheck found: ${CPPCHECK}")
  list(
    APPEND CPPCHECK_ARGS
    "--enable=all"
    "--inconclusive"
    "--quiet"
    "--std=c11"
    "--suppress=unusedFunction"
    "--suppress=missingInclude"
  )
  add_custom_target(
    cppcheck
    COMMAND ${CPPCHECK} ${CPPCHECK_ARGS} ${BUILD_FILES}
  )
else()
  message(STATUS "cppcheck not found")
endif()

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  message(STATUS "clang-format found: ${CLANG_FORMAT}")
  add_custom_target(
    clang-format
    COMMAND ${CLANG_FORMAT} -i --verbose -style=Google ${BUILD_FILES}
  )
else()
  message(STATUS "clang-format not found")
endif()