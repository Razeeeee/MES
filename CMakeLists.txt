cmake_minimum_required(VERSION 3.16)
project(MES_Project)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP
find_package(OpenMP REQUIRED)

# Set build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fopenmp")

# Include directories
include_directories(src)

# Source files - explicitly list to avoid archive folder
file(GLOB SOURCES 
    "src/*.cpp" 
    "src/*.h"
)

# Create executable
add_executable(MES_Project ${SOURCES})

# Link OpenMP
target_link_libraries(MES_Project PRIVATE OpenMP::OpenMP_CXX)

# Set output directory
set_target_properties(MES_Project PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
