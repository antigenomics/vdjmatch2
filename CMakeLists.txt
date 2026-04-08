cmake_minimum_required(VERSION 3.16)
project(vdjmatch2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_executable(vdjmatch2
        main/main.cpp
        src/CliConfig.cpp
        src/TsvReader.cpp
        src/Trie.cpp
        src/MatchWriter.cpp
        src/RepertoireMatcher.cpp
)

target_compile_options(vdjmatch2 PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/Ox /fp:fast>
        $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3 -ffast-math>
)

target_include_directories(vdjmatch2 PRIVATE include)
target_link_libraries(vdjmatch2 PRIVATE Threads::Threads)

install(TARGETS vdjmatch2
        RUNTIME DESTINATION vdjmatch2
        COMPONENT python
)