cmake_minimum_required(VERSION 3.16)
project(MES_Project)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Include directories
include_directories(src)

# Source files - explicitly list all source files
set(SOURCES
    src/main.cpp
    src/Config.cpp
    src/Element.cpp
    src/EquationSystem.cpp
    src/GaussianIntegration.cpp
    src/GlobalData.cpp
    src/Grid.cpp
    src/GridSelector.cpp
    src/Node.cpp
    src/Utils.cpp
)

# Create executable
add_executable(MES_Project ${SOURCES})

# Set output directory
set_target_properties(MES_Project PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
