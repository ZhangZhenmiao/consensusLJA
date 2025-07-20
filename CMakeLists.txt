cmake_minimum_required(VERSION 3.16)
project(cLJA VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED ENV{CONDA_PREFIX})
    message(FATAL_ERROR "Conda environment not activated! Please run:\n  conda activate pangaea")
endif()

message(STATUS "Using Conda environment at: $ENV{CONDA_PREFIX}")

# === Set Conda-specific paths ===
set(CONDA_PREFIX "$ENV{CONDA_PREFIX}")
set(CONDA_LIB_DIR "${CONDA_PREFIX}/lib")
set(CONDA_INCLUDE_DIR "${CONDA_PREFIX}/include")

# Make sure CMake searches Conda directories
list(APPEND CMAKE_PREFIX_PATH "${CONDA_PREFIX}")
list(APPEND CMAKE_LIBRARY_PATH "${CONDA_LIB_DIR}")
list(APPEND CMAKE_INCLUDE_PATH "${CONDA_INCLUDE_DIR}")

# Make linker prefer Conda libs
link_directories("${CONDA_LIB_DIR}")

# === Find required packages ===
find_package(OpenMP REQUIRED)

find_library(HTS_LIB hts PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(ZLIB_LIB z PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)
find_library(BZ2_LIB bz2 PATHS "${CONDA_LIB_DIR}" NO_DEFAULT_PATH REQUIRED)

find_path(HTS_INCLUDE_DIR hts.h PATHS "${CONDA_INCLUDE_DIR}/htslib" NO_DEFAULT_PATH REQUIRED)
find_path(ZLIB_INCLUDE_DIR zlib.h PATHS "${CONDA_INCLUDE_DIR}" NO_DEFAULT_PATH REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${ZLIB_INCLUDE_DIR}")

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

# Standalone tool for chimeric analysis
add_executable(analyze_chimeric src/scripts/analyze_chimeric.cpp)

# Include directories
target_include_directories(cLJA PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmdline
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/edlib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src # Assuming you might need this
)

# Link libraries
target_link_libraries(cLJA PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(analyze_chimeric PRIVATE ${HTS_LIB} ${ZLIB_LIB} ${BZ2_LIB})

# Compile options
target_compile_features(cLJA PRIVATE cxx_std_20)
target_compile_options(cLJA PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g>
    ${OpenMP_CXX_FLAGS}
)

# Better to explicitly link OpenMP flags
target_link_options(cLJA PRIVATE ${OpenMP_CXX_FLAGS})