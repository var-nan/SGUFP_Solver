cmake_minimum_required(VERSION 3.27)
project(SGUFP_Solver)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)


add_library(gurobi_c++ STATIC IMPORTED)
set_property(TARGET gurobi_c++ PROPERTY IMPORTED_LOCATION "/opt/gurobi1102/linux64/lib/libgurobi_c++.a")

add_library(gurobi110 SHARED IMPORTED)
set_property(TARGET gurobi110 PROPERTY IMPORTED_LOCATION "/opt/gurobi1102/linux64/lib/libgurobi110.so")

add_executable(SGUFP_Solver main.cpp
        DD.h
        Network.h
        Network.cpp
        DD.cpp
        Solver.cpp
        Solver.h
        grb.cpp
        grb.h)

# target_compile_features(SGUFP_Solver PUBLIC cxx_std_17)
target_link_libraries(SGUFP_Solver PRIVATE gurobi_c++)
target_link_libraries(SGUFP_Solver PRIVATE gurobi110)
target_include_directories(SGUFP_Solver PRIVATE "/opt/gurobi1102/linux64/include/")