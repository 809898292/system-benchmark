cmake_minimum_required(VERSION 3.10)
project(SystemBenchmark)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required libraries
find_package(ZLIB REQUIRED)
find_package(OpenSSL REQUIRED)

# Add executable
add_executable(SystemBenchmark main.cpp)

# Link libraries
target_link_libraries(SystemBenchmark 
    PRIVATE 
    ZLIB::ZLIB 
    OpenSSL::SSL 
    OpenSSL::Crypto
    m
)

# Enable optimizations
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(SystemBenchmark PRIVATE -O3 -march=native)
endif()

# Installation
install(TARGETS SystemBenchmark DESTINATION bin)