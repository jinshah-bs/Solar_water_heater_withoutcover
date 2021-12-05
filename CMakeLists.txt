cmake_minimum_required(VERSION 3.20)
project(Solar_analysis)

set(CMAKE_CXX_STANDARD 14)

add_executable(Solar_analysis main.cpp read_data.cpp read_data.h mat_prop.cpp mat_prop.h TDMA.cpp TDMA.h HTC.cpp HTC.h air_properties.h Write_data.cpp Write_data.h Solar_data.cpp Solar_data.h)
