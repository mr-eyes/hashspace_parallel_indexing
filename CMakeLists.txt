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

include_directories("${PROJECT_SOURCE_DIR}/lib/progressbar/include")

# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------

add_executable (index modified_index.cpp index.cpp lib.cpp)
target_link_libraries (index kProcessor z)
target_include_directories(index INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (pairwise index_and_pairwise.cpp modified_index.cpp kSpider_pairwise.cpp lib.cpp)
target_link_libraries (pairwise kProcessor z)
target_include_directories(pairwise INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (pairwise_inmemory lib.cpp inmemory_index_and_pairwise.cpp modified_index.cpp kSpider_pairwise.cpp)
target_link_libraries (pairwise_inmemory kProcessor z)
target_include_directories(pairwise_inmemory INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (index_and_pairwise_single lib.cpp index_and_pairwise_single.cpp modified_index.cpp kSpider_pairwise.cpp)
target_link_libraries (index_and_pairwise_single kProcessor z)
target_include_directories(index_and_pairwise_single INTERFACE ${kProcessor_INCLUDE_PATH})


add_executable (ranged_index lib.cpp ranged_index.cpp)
target_link_libraries (ranged_index kProcessor z)
target_include_directories(ranged_index INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (histo histogram.cpp lib.cpp)
target_link_libraries (histo kProcessor z)
target_include_directories(histo INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (bin_compare bin_pairwise.cpp lib.cpp)
target_link_libraries (bin_compare kProcessor z)
target_include_directories(bin_compare INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (ranged_pairwise lib.cpp ranged_pairwise.cpp)
target_link_libraries (ranged_pairwise kProcessor z)
target_include_directories(ranged_pairwise INTERFACE ${kProcessor_INCLUDE_PATH})