cmake_minimum_required (VERSION 3.4)
project (kTests C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC -lgomp -lrt -fopenmp -O3")


# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

set(kProcessor "${PROJECT_SOURCE_DIR}/lib/kProcessor")
add_subdirectory(${kProcessor})
set(kProcessor_INCLUDE_PATH "${kProcessor}/include")

# Add include directories.
include_directories("${kProcessor}/ThirdParty/MQF/include")
include_directories("${kProcessor}/ThirdParty/Blight")
include_directories (${kProcessor_INCLUDE_PATH})
include_directories("${kProcessor}/ThirdParty/caches/include")


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

add_executable (index index.cpp)
target_link_libraries (index kProcessor z)
target_include_directories(index INTERFACE ${kProcessor_INCLUDE_PATH})