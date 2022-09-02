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
include_directories("${PROJECT_SOURCE_DIR}/lib/json_parser/lib/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/zstr/src")



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


add_executable (ranged_index lib.cpp modified_index.cpp ranged_index.cpp)
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

add_executable (combine_pairwise lib.cpp combine_pairwise.cpp)
target_link_libraries (combine_pairwise kProcessor z)
target_include_directories(combine_pairwise INTERFACE ${kProcessor_INCLUDE_PATH})



add_executable (inmemory_batch_splitted_index lib.cpp modified_index.cpp inmemory_batch_splitted_index.cpp)
target_link_libraries (inmemory_batch_splitted_index kProcessor z)
target_include_directories(inmemory_batch_splitted_index INTERFACE ${kProcessor_INCLUDE_PATH})


add_executable (sig_to_bin sig_to_bin.cpp)
target_link_libraries (sig_to_bin kProcessor z)
target_include_directories(sig_to_bin INTERFACE ${kProcessor_INCLUDE_PATH})

add_executable (sigs_to_bins sigs_to_bins.cpp)
target_link_libraries (sigs_to_bins kProcessor z)
target_include_directories(sigs_to_bins INTERFACE ${kProcessor_INCLUDE_PATH})