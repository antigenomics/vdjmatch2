cmake_minimum_required(VERSION 3.16)
project(vdjmatch2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include)

add_executable(vdjmatch2
        main/main.cpp
        src/CliConfig.cpp
        src/TsvReader.cpp
        src/Trie.cpp
        src/MatchWriter.cpp
        src/RepertoireMatcher.cpp
)


enable_testing()

add_executable(vdjmatch2_tests
        tests/run_tests.cpp
        src/CliConfig.cpp
        src/TsvReader.cpp
        src/Trie.cpp
        src/MatchWriter.cpp
        src/RepertoireMatcher.cpp
)

add_test(NAME vdjmatch2_tests COMMAND vdjmatch2_tests)
