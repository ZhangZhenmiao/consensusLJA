cmake_minimum_required(VERSION 3.16)
project(cLJA VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Debug if no build type is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Create output directories before they're needed
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) # Changed from SOURCE_DIR to BINARY_DIR

# Find OpenMP
find_package(OpenMP REQUIRED)

# Collect source files properly
file(GLOB DBG_SOURCES src/dbg/*.cpp)
file(GLOB MULTIDBG_SOURCES src/multidbg/*.cpp)

# Add the executable
add_executable(cLJA
    src/main.cpp
    src/utils.cpp
    ${DBG_SOURCES}
    ${MULTIDBG_SOURCES}
    lib/edlib/src/edlib.cpp
)

# Include directories
target_include_directories(cLJA PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmdline
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/edlib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src # Assuming you might need this
)

# Link libraries
target_link_libraries(cLJA PRIVATE OpenMP::OpenMP_CXX)

# Compile options
target_compile_features(cLJA PRIVATE cxx_std_20)
target_compile_options(cLJA PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g>
    ${OpenMP_CXX_FLAGS}
)

# Better to explicitly link OpenMP flags
target_link_options(cLJA PRIVATE ${OpenMP_CXX_FLAGS})