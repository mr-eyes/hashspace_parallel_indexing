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
include_directories("${PROJECT_SOURCE_DIR}/include")

# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

add_executable (index modified_index.cpp index.cpp)
target_link_libraries (index kProcessor z)
target_include_directories(index INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (pairwise index_and_pairwise.cpp modified_index.cpp kSpider_pairwise.cpp)
target_link_libraries (pairwise kProcessor z)
target_include_directories(pairwise INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (histo histogram.cpp)
target_link_libraries (histo kProcessor z)
target_include_directories(histo INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (bin_compare bin_pairwise.cpp)
target_link_libraries (bin_compare kProcessor z)
target_include_directories(bin_compare INTERFACE ${kProcessor_INCLUDE_PATH})