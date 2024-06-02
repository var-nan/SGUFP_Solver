cmake_minimum_required(VERSION 3.27)
project(SGUFP_Solver)

set(CMAKE_CXX_STANDARD 17)

add_executable(SGUFP_Solver main.cpp
        Node.h
        DD.h
        Network.h
        Network.cpp)
