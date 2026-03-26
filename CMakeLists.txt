cmake_minimum_required(VERSION 3.16)
project(vdjmatch2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
include(CTest)

add_executable(vdjmatch2
        main/main.cpp
        src/CliConfig.cpp
        src/TsvReader.cpp
        src/Trie.cpp
        src/MatchWriter.cpp
        src/RepertoireMatcher.cpp
)

target_include_directories(vdjmatch2 PRIVATE include)

install(TARGETS vdjmatch2
        RUNTIME DESTINATION vdjmatch2
        COMPONENT python
)

if(BUILD_TESTING)
    add_executable(vdjmatch2_tests
            tests/run_tests.cpp
            src/CliConfig.cpp
            src/TsvReader.cpp
            src/Trie.cpp
            src/MatchWriter.cpp
            src/RepertoireMatcher.cpp
    )

    target_include_directories(vdjmatch2_tests PRIVATE include)

    add_test(NAME vdjmatch2_tests COMMAND vdjmatch2_tests)
endif()
