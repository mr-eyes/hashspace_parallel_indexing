#include "algorithms.hpp"
#include <parallel_hashmap/phmap.h>
#include <lib.hpp>


int main(int argc, char** argv) {

    if (argc < 4) {
        cout << "./index <bins_dir> <out_dir> <min_scale> <cores>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string out_dir = argv[2];
    int min_scale = stoi(argv[3]);
    int cores = stoi(argv[4]);
    uint64_t max_hash = UINT64_MAX / (uint64_t)min_scale;

    auto ranges = splitted_ranges(max_hash, cores);
    flat_hash_map<int, colored_kDataFrame*> frames_map;
    frames_map.reserve(cores);
    int kSize = 21;

    for (int i = 0; i < cores; i++) {
        frames_map.insert(pair(i, new colored_kDataFrame()));
    }


#pragma omp parallel num_threads(cores)
    {
#pragma omp for

        for (int i = 0; i < cores; i++) {
            uint64_t from_hash = get<0>(ranges[i]);
            uint64_t to_hash = get<1>(ranges[i]);

            cout << "INDEXING PART " << i + 1 << endl;
            bins_indexing_hashsplit(frames_map[i], bins_dir, kSize, from_hash, to_hash);
            frames_map[i]->save(out_dir + "/" + to_string(i) + "_");

            // tmp dump kmers to understand how it will works
            ofstream phmapDump(out_dir + "/" + to_string(i) + "_phmap_dump");
            auto it = frames_map[i]->getkDataFrame()->begin();
            while (it != frames_map[i]->getkDataFrame()->end())
            {
                phmapDump << it.getHashedKmer() << "," << it.getCount() << endl;
                it++;
            }
            phmapDump.close();

            // frames_map.erase(i);
        }
    }
}

