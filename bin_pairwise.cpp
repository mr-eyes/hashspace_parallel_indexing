#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include <ctime>
#include<omp.h>
#include <glob.h>
#include <string>
#include <stdexcept>
#include "parallel_hashmap/phmap_dump.h"
#include <cstdlib>

using namespace std;
// using namespace phmap;


int main(int argc, char** argv) {

    if (argc != 3) {
        cout << "run: ./check_bin <bin1> <bin2>" << endl;
        exit(1);
    }

    string bin1_path = argv[1];
    string bin2_path = argv[2];

    phmap::flat_hash_set<uint64_t> bin1_hashes;
    phmap::BinaryInputArchive ar1_in(bin1_path.c_str());
    bin1_hashes.phmap_load(ar1_in);

    phmap::flat_hash_set<uint64_t> bin2_hashes;
    phmap::BinaryInputArchive ar2_in(bin2_path.c_str());
    bin2_hashes.phmap_load(ar2_in);

    uint64_t shared_hashes = count_if(bin1_hashes.begin(), bin1_hashes.end(), [&](uint64_t k) {return bin2_hashes.find(k) != bin2_hashes.end();});

    cout << "size(" << bin1_path << "):  " << bin1_hashes.size() << endl;
    cout << "size(" << bin2_path << "):  " << bin2_hashes.size() << endl; 
    cout << "shared_kmers: " << shared_hashes << endl;
    float max_containment = (float)shared_hashes / min(bin1_hashes.size(), bin2_hashes.size());
    cout << "containment: " << to_string(max_containment) << endl;



}