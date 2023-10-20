cmake_minimum_required(VERSION 3.13)
project(MathUtils)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.h
)

set(SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.c
)

add_library(${PROJECT_NAME} STATIC ${HEADERS} ${SOURCES})

set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)

target_compile_options(
  ${PROJECT_NAME}
  PRIVATE
  -Wall
  -Werror
  -Wextra
  -Wpedantic
  -std=c11
  -D_GNU_SOURCE
)

target_link_libraries(${PROJECT_NAME} PRIVATE m)

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
    COMMAND ${CPPCHECK} ${CPPCHECK_ARGS} ${HEADERS} ${SOURCES}
  )
else()
  message(STATUS "cppcheck not found")
endif()

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  message(STATUS "clang-format found: ${CLANG_FORMAT}")
  add_custom_target(
    clang-format
    COMMAND ${CLANG_FORMAT} -i --verbose -style=Google ${HEADERS} ${SOURCES}
  )
else()
  message(STATUS "clang-format not found")
endif()