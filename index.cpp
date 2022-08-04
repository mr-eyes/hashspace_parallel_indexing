#include "algorithms.hpp"
#include <parallel_hashmap/phmap.h>

void bins_indexing_hashsplit(colored_kDataFrame* res, string bins_dir, int selective_kSize, uint64_t from_hash, uint64_t to_hash);

vector<tuple<uint64_t, uint64_t>> splitted_ranges(uint64_t max_hash, int cores) {
    vector<tuple<uint64_t, uint64_t>> ranges;
    uint64_t from_hash = 0;
    uint64_t to_hash = 0;
    uint64_t step = (uint64_t)(max_hash / cores);
    for (int i = 0; i < cores; i++) {
        to_hash += step;
        ranges.push_back({ from_hash, to_hash });
        from_hash += step;
    }
    ranges[cores - 1] = { from_hash, max_hash };
    return ranges;
}

int main(int argc, char** argv) {

    if(argc < 4){
        cout << "./index <bins_dir> <out_dir> <cores>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string out_dir = argv[2];
    int cores = stoi(argv[3]);

    auto ranges = splitted_ranges(18446744073709551, cores);
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

            auto it = frames_map[i]->frame->begin();
            while(it != frames_map[i]->frame->end())
            {
                phmapDump << it.getHashedKmer() << ","<< it.getCount()<<endl;
                it++;
            }
            phmapDump.close();

            // frames_map.erase(i);
        }
    }
}

