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
#include <ctime>
#include "progressbar.hpp"

#ifndef LOGGING 
#define LOGGING 0
#endif 

#include <lib.hpp>

using namespace std;

int main(int argc, char** argv) {

    if (argc < 6) {
        cout << "run: ./inmemory_batch_splitted_index <bins_dir> <output_prefix> <total_parts> <scale> <cores>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string output_prefix = argv[2];
    uint64_t total_parts = to_uint64_t(argv[3]);
    uint64_t scale = to_uint64_t(argv[4]);
    int cores = stoi(argv[5]);

    if (cores > total_parts) {
        cout << "downscaling cores to " << total_parts << endl;
        cores = total_parts;
    }

    /*
    multi-threaded loading of all signatures and store them in memory.
    */

    // Load all files in memory
    auto* bin_to_hashes = new BINS_PHMAP();
    auto* groupName_to_kmerCount = new BINS_KMER_COUNT();
    flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
    load_all_bins(bins_dir, bin_to_hashes, groupName_to_kmerCount, cores);

    uint64_t max_hash = UINT64_MAX / (uint64_t)scale;
    auto ranges = splitted_ranges(max_hash, total_parts);

    auto begin_time = Time::now();

    progressbar totalbar(cores);
    totalbar.set_todo_char(" ");
    totalbar.set_done_char("█");
    totalbar.set_opening_bracket_char("[");
    totalbar.set_closing_bracket_char("]");

    int thread_num_1, num_threads_1, start_1, end_1, vec_i_1;
    int n_1 = ranges.size();
    omp_set_num_threads(cores);
    bool _names_to_id_map_flag_ = false;


#pragma omp parallel private(vec_i_1,thread_num_1,num_threads_1,start_1,end_1)
    {

        thread_num_1 = omp_get_thread_num();
        num_threads_1 = omp_get_num_threads();
        start_1 = thread_num_1 * n_1 / num_threads_1;
        end_1 = (thread_num_1 + 1) * n_1 / num_threads_1;

        for (vec_i_1 = start_1; vec_i_1 != end_1; ++vec_i_1) {
            uint64_t from_hash = get<0>(ranges[vec_i_1]);
            uint64_t to_hash = get<1>(ranges[vec_i_1]);
            flat_hash_map<string, uint32_t> name_to_id_map;

            ranged_index_new(bin_to_hashes,
                output_prefix,
                from_hash,
                to_hash,
                vec_i_1 + 1,
                name_to_id_map);
#pragma omp critial
            {
                if (_names_to_id_map_flag_ == false) {
                    for (auto& [_NAME, _COUNT] : *groupName_to_kmerCount) {
                        groupID_to_kmerCount[name_to_id_map[_NAME]] = _COUNT;
                    }

                    std::ofstream fstream_id_to_name;
                    fstream_id_to_name.open(output_prefix + "_id_to_name.tsv");
                    fstream_id_to_name << "ID\tname\n";
                    for (auto& [_name, _id] : name_to_id_map) {
                        fstream_id_to_name << _id << '\t' << _name << '\n';
                    }
                    fstream_id_to_name.close();

                    phmap::BinaryOutputArchive ar_out(string(output_prefix + "_groupID_to_kmerCount.bin").c_str());
                    groupID_to_kmerCount.phmap_dump(ar_out);

                    delete groupName_to_kmerCount;
                    _names_to_id_map_flag_ = true;

                }

            }
        }
#pragma omp critial
        totalbar.update();
    }

    /*
    END loading.
    */

    cout << endl;

    std::ofstream metadata;
    metadata.open(output_prefix + ".metadata");
    metadata << "scale," << to_string(scale) << "\n";
    metadata << "parts," << to_string(total_parts) << "\n";
    metadata.close();

#ifdef LOGGING
    auto total = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
    cout << "======================\n DONE PROCESSING " << total_parts << " parts in " << total << " secs." << endl;
#endif
}