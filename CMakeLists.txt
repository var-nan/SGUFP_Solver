cmake_minimum_required(VERSION 3.27)
project(SGUFP_Solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(SOURCES Network.cpp DD.cpp grb.cpp OriginalProblem.h)
add_executable(SGUFP_Solver main.cpp ${SOURCES})

# preprocessor macros.
add_compile_definitions(MAX_WIDTH=5)

# google test dependency
include(FetchContent)
FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()


if (UNIX)

    message(STATUS "Current Build type: ${CMAKE_BUILD_TYPE}")

    # gurobi libraries
    set(GUROBI_HOME "/opt/gurobi1102/linux64/")
    set(GUROBI_LIB "${GUROBI_HOME}lib/")
    set(GUROBI_INCLUDE "${GUROBI_HOME}include/")

    add_library(gurobi_c++ STATIC IMPORTED)
    set_property(TARGET gurobi_c++ PROPERTY IMPORTED_LOCATION "${GUROBI_LIB}libgurobi_c++.a")

    add_library(gurobi110 SHARED IMPORTED)
    set_property(TARGET gurobi110 PROPERTY IMPORTED_LOCATION "${GUROBI_LIB}libgurobi110.so")

    add_compile_definitions(DEBUG)
    target_compile_options(SGUFP_Solver PRIVATE -Wall -Wextra -O0 ) # add -O2 flag in final phase. -Wunused
    # compiler flags: -Wuseless-cast -Wshadow, -Wduplicated-cond, -Wduplicated-branches
    # compiler flags: -fsanitize=pointer-subtract, -fsanitize=address
    # compiler flags : -Wnull-dereference, -Wold-style-cast
    # add -Werror flag later
    target_compile_options(SGUFP_Solver PRIVATE -Wuseless-cast -Wnull-dereference)
    # target_compile_features(SGUFP_Solver PUBLIC cxx_std_17)
    target_link_libraries(SGUFP_Solver PRIVATE gurobi_c++)
    target_link_libraries(SGUFP_Solver PRIVATE gurobi110)
    target_include_directories(SGUFP_Solver PRIVATE ${GUROBI_INCLUDE})

    # test targets
    add_executable(tests2 tests/tests.cpp ${SOURCES})
    target_link_libraries(tests2 gurobi_c++ gurobi110 GTest::gtest_main)
    target_include_directories(tests2 PRIVATE ${GUROBI_INCLUDE})

    add_executable(testmain testmain.cpp ${SOURCES})
    target_link_libraries(testmain gurobi_c++ gurobi110 GTest::gtest_main)
    target_include_directories(testmain PRIVATE ${GUROBI_INCLUDE})

    include(GoogleTest)
    gtest_discover_tests(testmain)
    gtest_discover_tests(tests2)

elseif(WIN32)
    # windows system
    option(CXX "enable C++ compilation" ON)

    if (CXX)
        enable_language(CXX)
    endif()

    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
    if (CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Debug")
    endif()

    if (MSVC AND MT)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
    endif()

    find_package(GUROBI REQUIRED)

    include_directories(${GUROBI_INCLUDE_DIRS})

    add_compile_definitions(DEBUG)

    message(STATUS "GUROBI_CXX_LIBRARY: ${GUROBI_CXX_LIBRARY}")
    message(STATUS "GUROBI_CXX_DEBUG LIBRARY: ${GUROBI_CXX_DEBUG_LIBRARY}")

    target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
            debug ${GUROBI_CXX_DEBUG_LIBRARY})

    target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})

    if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
        include(FeatureSummary)
        feature_summary(WHAT ALL)
    endif()

    add_executable(testmain testmain.cpp ${SOURCES})
    target_link_libraries(testmain GTest::gtest_main)

    add_executable(tests2 tests/tests.cpp ${SOURCES})
    target_link_libraries(tests2 GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(testmain)
    gtest_discover_tests(tests2)

endif()