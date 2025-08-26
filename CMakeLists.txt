cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0074 NEW)
project(cLJA VERSION 1.0 LANGUAGES CXX)

# === C++ standard ===
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# === Dependencies ===
find_package(OpenMP REQUIRED)
if(DEFINED ENV{CONDA_PREFIX})
    set(ZLIB_ROOT $ENV{CONDA_PREFIX})
endif()
find_package(ZLIB REQUIRED)

include_directories(${ZLIB_INCLUDE_DIRS})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# === Output directories ===
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Ensure directories exist
file(MAKE_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

# === Sources ===
file(GLOB DBG_SOURCES src/dbg/*.cpp)
file(GLOB MULTIDBG_SOURCES src/multidbg/*.cpp)

add_executable(cLJA
    src/main.cpp
    src/utils.cpp
    ${DBG_SOURCES}
    ${MULTIDBG_SOURCES}
    lib/edlib/src/edlib.cpp
)

# === Includes ===
target_include_directories(cLJA PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/cmdline
    ${CMAKE_CURRENT_SOURCE_DIR}/lib/edlib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# === Libraries ===
target_link_libraries(cLJA PRIVATE OpenMP::OpenMP_CXX)

# === Compiler flags ===
target_compile_features(cLJA PRIVATE cxx_std_20)
target_compile_options(cLJA PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-g>
)

# === Linker flags (OpenMP explicitly) ===
target_link_options(cLJA PRIVATE ${OpenMP_CXX_FLAGS})

# === Submodules ===
add_subdirectory(lib/LJA)
