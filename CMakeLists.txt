cmake_minimum_required(VERSION 3.10)
project(MyProject)

# Add Eigen directory
include_directories(${CMAKE_SOURCE_DIR}/eigen)

# Add ThreadManager directory
include_directories(${CMAKE_SOURCE_DIR}/ThreadManager)

# Add your source files
add_executable(MyProject extension.cpp matrix.cpp ThreadManager/ThreadManager.cpp)

# Add libraries
target_link_libraries(MyProject pthread)
