cmake_minimum_required(VERSION 3.7)
project(EKF_Project)

set(CMAKE_CXX_STANDARD 11)

set(sources
        src/FusionEKF.cpp
        src/kalman_filter.cpp
        src/main.cpp
        src/tools.cpp)

add_executable(EKF_Project ${sources})