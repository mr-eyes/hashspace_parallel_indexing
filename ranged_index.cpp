#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include "algorithms.hpp"
#include <vector>
#include <stdexcept>
#include <sstream>
#include <string>
#include <fstream>
#include "parallel_hashmap/phmap_dump.h"
#include <queue>
#include <lib.hpp>
#include <ctime>
#include "progressbar.hpp"

using namespace std;

int main(int argc, char** argv) {

    if (argc < 6) {
        cout << "run: ./ranged_index <bins_dir> <output_prefix> <this_part> <total_parts> <scale>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string output_prefix = argv[2];
    uint64_t this_part = to_uint64_t(argv[3]);
    uint64_t total_parts = to_uint64_t(argv[4]);
    uint64_t scale = to_uint64_t(argv[5]);
    int loading_cores = 16;

    // Load all files in memory
    auto* bin_to_hashes = new BINS_PHMAP();
    auto* groupName_to_kmerCount = new BINS_KMER_COUNT();
    load_all_bins(bins_dir, bin_to_hashes, groupName_to_kmerCount, loading_cores);


    if (this_part <= 0) {
        cout << "this_part must be >= 1" << endl;
        exit(1);
    }

    uint64_t max_hash = UINT64_MAX / (uint64_t)scale;
    auto ranges = splitted_ranges(max_hash, total_parts);
    uint64_t from_hash = get<0>(ranges[this_part - 1]);
    uint64_t to_hash = get<1>(ranges[this_part - 1]);

#if LOGGING == 1
    cout << "This run: " << endl;
    cout << "\t\t- from hash: " << to_string(from_hash) << endl;
    cout << "\t\t- to hash: " << to_string(to_hash) << endl;
    cout << "--------------------------" << endl;
#endif

    auto begin_time = Time::now();
    flat_hash_map<string, uint32_t> groupNameMap;
    ranged_index_new(bin_to_hashes, output_prefix, from_hash, to_hash, this_part, groupNameMap);
#if LOGGING == 1
    auto total = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
    cout << "======================\n DONE PROCESSING PART " << this_part << "/" << total_parts << " in " << total << " secs." << endl;
#endif
}