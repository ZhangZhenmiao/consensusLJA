project(sequences)
set(CMAKE_CXX_STANDARD 17)

include_directories(.)
add_library(lja_sequence STATIC contigs.cpp sequence.cpp seqio.cpp)

target_link_libraries (lja_sequence ${CMAKE_THREAD_LIBS_INIT} ${ZLIB_LIBRARIES} m)
target_include_directories (lja_sequence PUBLIC ${ZLIB_INCLUDE_DIRS})