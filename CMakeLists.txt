# FILEPATH: /mnt/e/Projects/emc/CMakeLists.txt

cmake_minimum_required(VERSION 3.0)
project(emc)

# Include Torch
list(APPEND CMAKE_PREFIX_PATH "/home/w/libtorch")
find_package(Torch REQUIRED)
include_directories(${TORCH_INCLUDE_DIRS})

# Add your project files and targets here
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")
add_executable(emc mc.cpp)

target_link_libraries(emc "${TORCH_LIBRARIES}")
set_property(TARGET emc PROPERTY CXX_STANDARD 17)
