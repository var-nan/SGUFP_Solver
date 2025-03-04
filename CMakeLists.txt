cmake_minimum_required(VERSION 3.27)
project(SGUFP_Solver)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)

set(SOURCES Network.cpp DD.cpp Cut.h grb.cpp DDSolver.cpp NodeExplorer.cpp optimized.cpp OriginalProblem.h)
#set(SOURCES Network.cpp DD.cpp NodeExplorer.cpp grb.cpp)
add_executable(SGUFP_Solver main.cpp ${SOURCES})

# preprocessor macros.
add_compile_definitions(MAX_WIDTH=32)


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Fetching google test")
    # google test dependency
    include(FetchContent)
    FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    enable_testing()
endif()

if (UNIX)
    message(STATUS "Current Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "CMake Version: ${CMAKE_VERSION}")
    message(STATUS "CMake generator: ${CMAKE_GENERATOR}")
    message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
    message(STATUS "C compiler version: ${CMAKE_C_COMPILER_VERSION}")
    message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
    message(STATUS "C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")

    # set build specific compiler flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funroll-loops -fno-rtti -fno-exceptions")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -fno-omit-frame-pointer -fno-rtti")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")

    # gurobi libraries
    set(GUROBI_HOME "/opt/gurobi1103/linux64/")
    set(GUROBI_LIB "${GUROBI_HOME}lib/")
    set(GUROBI_INCLUDE "${GUROBI_HOME}include/")

    add_library(gurobi_c++ STATIC IMPORTED)
    set_property(TARGET gurobi_c++ PROPERTY IMPORTED_LOCATION "${GUROBI_LIB}libgurobi_c++.a")

    add_library(gurobi110 SHARED IMPORTED)
    set_property(TARGET gurobi110 PROPERTY IMPORTED_LOCATION "${GUROBI_LIB}libgurobi110.so")

    add_compile_definitions(SOLVER_STATS)
    target_compile_options(SGUFP_Solver PRIVATE -Wall -Wextra -pg -g -pthread) # add -O3 flag in final phase. -Wunused
    target_link_options(SGUFP_Solver PRIVATE -pg -g -pthread)

    find_package(GSL REQUIRED)
    target_link_libraries(SGUFP_Solver PRIVATE GSL::gsl)
#    target_compile_options(SGUFP_Solver PRIVATE -Wall -Wextra -Wpedantic -Werror -Wshadow -Wnull-dereference)
#    target_compile_options(SGUFP_Solver PRIVATE -Wuseless-cast -Wduplicated-cond -Wduplicated-branches -Wold-style-cast)
#    target_compile_options(SGUFP_Solver PRIVATE -fsanitize=pointer-subtract -fsanitize=address)
    # compiler flags: -Wuseless-cast -Wshadow, -Wduplicated-cond, -Wduplicated-branches
    # compiler flags: -fsanitize=pointer-subtract, -fsanitize=address -fno-rtti,
    # compiler flags : -Wnull-dereference, -Wold-style-cast
    # add -Werror flag later
#    target_compile_options(SGUFP_Solver PRIVATE -Wuseless-cast -Wnull-dereference)
    # target_compile_features(SGUFP_Solver PUBLIC cxx_std_17)
    target_link_libraries(SGUFP_Solver PRIVATE gurobi_c++)
    target_link_libraries(SGUFP_Solver PRIVATE gurobi110)
    target_include_directories(SGUFP_Solver PRIVATE ${GUROBI_INCLUDE})

#    message(STATUS "General FLAGS: ${CMAKE_CXX_FLAGS}")
    message(STATUS "RELEASE FLAGS: ${CMAKE_CXX_FLAGS_RELEASE}")
    message(STATUS "DEBUG FLAGS: ${CMAKE_CXX_FLAGS_DEBUG}")
    message(STATUS "RELWITHDEBINFO FLAGS: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")


    # test targets
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_executable(tests2 tests/tests.cpp ${SOURCES}
                tests/test_utils.h)
        target_link_libraries(tests2 gurobi_c++ gurobi110 GTest::gtest_main)
        target_include_directories(tests2 PRIVATE ${GUROBI_INCLUDE})

        add_executable(testmain testmain.cpp ${SOURCES})
        target_link_libraries(testmain gurobi_c++ gurobi110 GTest::gtest_main)
        target_include_directories(testmain PRIVATE ${GUROBI_INCLUDE})


    #    if (CMAKE_BUILD_TYPE STREQUAL "Release")
    #        set_target_properties(testmain PROPERTIES EXCLUDE_FROM_ALL TRUE)
    #        set_target_properties(tests2 PROPERTIES EXCLUDE_FROM_ALL TRUE)
    #    endif()
        include(GoogleTest)
        gtest_discover_tests(testmain)
        gtest_discover_tests(tests2)
    endif()

elseif(WIN32)
    # windows system
    option(CXX "enable C++ compilation" ON)

    if (CXX)
        enable_language(CXX)
    endif()

    list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
#    if (CMAKE_BUILD_TYPE)
#        set(CMAKE_BUILD_TYPE "Debug")
#    endif()

    message(STATUS "Current Build type: ${CMAKE_BUILD_TYPE}")
    message(STATUS "CMake Version: ${CMAKE_VERSION}")
    message(STATUS "CMake generator: ${CMAKE_GENERATOR}")
    message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
    message(STATUS "C compiler version: ${CMAKE_C_COMPILER_VERSION}")
    message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
    message(STATUS "C++ compiler version: ${CMAKE_CXX_COMPILER_VERSION}")


#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp /openmp:llvm")
#
#    if (MSVC AND MT)
#        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT /openmp")
#        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd /openmp")
#        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT /openmp")
#        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd /openmp")
#    endif()
#    # flag specific for Dr Memory
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Zi /Oy- /MTd /openmp /openmp:llvm")
#    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /openmp /openmp:llvm")
##    message(STATUS "${CMAKE_CXX_FLAGS_DEBUG}")
    find_package(GUROBI REQUIRED)

    include_directories(${GUROBI_INCLUDE_DIRS})

#    add_compile_definitions(DEBUG)
    add_compile_definitions(SOLVER_STATS)

    find_package(GSL REQUIRED)
    target_link_libraries(SGUFP_Solver GSL::gsl)

#    target_compile_options(SGUFP_Solver PRIVATE -O0 -g -Wextra)
#    target_link_options(SGUFP_Solver PRIVATE -g -O0 -fsanitize=addr)

    message(STATUS "GUROBI_CXX_LIBRARY: ${GUROBI_CXX_LIBRARY}")
    message(STATUS "GUROBI_CXX_DEBUG LIBRARY: ${GUROBI_CXX_DEBUG_LIBRARY}")

    target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
            debug ${GUROBI_CXX_DEBUG_LIBRARY})

    target_link_libraries(${CMAKE_PROJECT_NAME} ${GUROBI_LIBRARY})

    if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
        include(FeatureSummary)
        feature_summary(WHAT ALL)
    endif()

#    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
#        add_executable(testmain testmain.cpp ${SOURCES})
#        target_link_libraries(testmain GTest::gtest_main)
#    message(STATUS "Fetching google test")
#    # google test dependency
#    include(FetchContent)
#    FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git
#            GIT_TAG v1.15.2
#    )
#    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#    FetchContent_MakeAvailable(googletest)
#    enable_testing()
#
#        add_executable(tests2 tests/tests2.cpp DD.cpp Network.cpp tests/test_utils.h)
#        target_link_libraries(tests2 optimized ${GUROBI_CXX_LIBRARY}
#                debug ${GUROBI_CXX_DEBUG_LIBRARY})
#        target_link_libraries(tests2 ${GUROBI_LIBRARY})
#        target_link_libraries(tests2 GTest::gtest_main)
#
#        include(GoogleTest)
##        gtest_discover_tests(testmain)
#        gtest_discover_tests(tests2)
##    endif()
endif()