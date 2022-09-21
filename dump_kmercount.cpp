#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <vector>
#include <iostream>
#include <lib.hpp>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;


int main(int argc, char** argv) {

    string input_prefix = argv[1];

    phmap::flat_hash_map<uint32_t, uint32_t> group_id_to_kmer_count;
    {
        phmap::BinaryInputArchive ar_in(string(input_prefix + "_groupID_to_kmerCount.bin").c_str());
        group_id_to_kmer_count.phmap_load(ar_in);
    }
    assert(group_id_to_kmer_count.size());

    std::ofstream result;
    result.open(input_prefix + "_groupID_to_kmerCount.tsv");

    for (auto& [k, v] : group_id_to_kmer_count) {
        result << k << "\t" << v << endl;
    }

    result.close();
}