cmake_minimum_required(VERSION 3.24)
include_guard(GLOBAL)

project(pyRxTools VERSION 1.1.0)

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

include(pyRxTools.cmake)

add_library(pyRxTools SHARED ${PYRXTOOLS_SOURCES})

add_custom_command(TARGET pyRxTools POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:pyRxTools>
        "${CMAKE_SOURCE_DIR}/../src/main/resources/base/bin/pyRxTools.dll"
)

target_include_directories(pyRxTools PRIVATE ${PYRXTOOLS_EXTERNAL_INCLUDES})
target_link_libraries(pyRxTools PRIVATE ${PYRXTOOLS_EXTERNAL_LINKS})
target_compile_definitions(pyRxTools PRIVATE PYRXTOOLS_EXPORTS)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT pyRxTools)

if (MSVC)
	target_compile_options(pyRxTools PRIVATE /W3 /WX)
else()
	target_compile_options(pyRxTools PRIVATE -Wall -Wextra -Werror)
endif()