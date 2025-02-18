cmake_minimum_required(VERSION 3.16)
project(cLJA VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Release if no build type is specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Add the executable
add_executable(cLJA
    src/main.cpp
    src/dot_graph.cpp
    src/read_aln.cpp
    src/graph_simplification.cpp
    lib/edlib/src/edlib.cpp
)

# Include directories
target_include_directories(cLJA PRIVATE
    lib/cmdline
    lib/edlib/include
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