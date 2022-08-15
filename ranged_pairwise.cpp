#include <iostream>
#include <cstdint>
#include <chrono>
#include "colored_kDataFrame.hpp"
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include <ctime>
#include "algorithms.hpp"
#include<omp.h>
#include <queue>
#include "parallel_hashmap/phmap_dump.h"
#include "lib.hpp"
#include "progressbar.hpp"
typedef std::chrono::high_resolution_clock Time;


void pairwise_ranged(
    PAIRS_COUNTER* edges,
    flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>* color_to_ids,
    flat_hash_map<uint64_t, uint64_t>& colorsCount,
    int user_threads
) {

    // convert map to vec for parallelization purposes.
    auto begin_col_to_ids = color_to_ids->begin();
    auto end_col_to_ids = color_to_ids->end();


    auto vec_color_to_ids = vector<pair<uint64_t, vector<uint32_t>>>();
    for (auto it : *color_to_ids) {
        vector<uint32_t> tmp(std::make_move_iterator(it.second.begin()), std::make_move_iterator(it.second.end()));
        vec_color_to_ids.push_back(make_pair(it.first, tmp));
    }

    int thread_num, num_threads, start, end, vec_i;
    int n = vec_color_to_ids.size();

    omp_set_num_threads(user_threads);

    progressbar bar((omp_get_thread_num() + 1) * n / omp_get_num_threads());
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");

#pragma omp parallel private(vec_i,thread_num,num_threads,start,end)
    {
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        start = thread_num * n / num_threads;
        end = (thread_num + 1) * n / num_threads;

        for (vec_i = start; vec_i != end; ++vec_i) {
            auto& item = vec_color_to_ids[vec_i];
            Combo combo = Combo();
            combo.combinations(item.second.size());
            for (uint64_t i = 0; i < combo.combs.size(); i++) {
                auto const& seq_pair = combo.combs[i];
                uint32_t _seq1 = item.second[seq_pair.first];
                uint32_t _seq2 = item.second[seq_pair.second];
                ascending(_seq1, _seq2);

                auto _p = make_pair(_seq1, _seq2);
                uint32_t ccount = colorsCount[item.first];


                edges->try_emplace_l(_p,
                    [ccount](PAIRS_COUNTER::value_type& v) { v.second += ccount; },           // called only when key was already present
                    ccount
                );
            }
#pragma omp critical
            bar.update();
        }
    }

}

void load_colors_to_sources(const std::string& filename, phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>* map)
{
    phmap::BinaryInputArchive ar_in(filename.c_str());
    size_t size;
    ar_in.loadBinary(&size);
    map->reserve(size);
    while (size--)
    {
        uint64_t k;
        phmap::flat_hash_set<uint32_t> v;
        ar_in.loadBinary(&k);
        ar_in.loadBinary(&v);
        map->insert_or_assign(std::move(k), std::move(v));
    }
}


int main(int argc, char** argv) {
    string input_prefix = argv[1];
    int part_id = stoi(argv[2]);
    int cores = stoi(argv[3]);

    string part_name = "part" + to_string(part_id);
    string _file_name_to_id = input_prefix + "_id_to_name.tsv";
    string _file_id_to_kmer_count = input_prefix + "_groupID_to_kmerCount.bin";
    string _file_color_to_sources = input_prefix + "_color_to_sources." + part_name + ".bin";
    string _file_color_to_count = input_prefix + "_color_count." + part_name + ".bin";

    // Loading all data
    cout << "Loading binaries ... ";
    auto begin_time = Time::now();
    auto* colors_to_sources = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>();
    load_colors_to_sources(_file_color_to_sources, colors_to_sources);
    assert(colors_to_sources->size());

    phmap::flat_hash_map<uint32_t, uint32_t> kmerCount;
    phmap::BinaryInputArchive ar_in_kmer_count(_file_id_to_kmer_count.c_str());
    kmerCount.phmap_load(ar_in_kmer_count);
    assert(kmerCount.size());

    flat_hash_map<uint64_t, uint64_t> colorsCount;
    phmap::BinaryInputArchive ar_in_colorsCount(_file_color_to_count.c_str());
    colorsCount.phmap_load(ar_in_colorsCount);
    assert(colorsCount.size());
    auto load_time = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
    cout << "Done in " << load_time << " secs." << endl << endl;
    // Loading done

    auto* edges = new PAIRS_COUNTER();

    begin_time = Time::now();
    cout << "performing pairwise comparisons..." << endl;
    pairwise_ranged(edges, colors_to_sources, colorsCount, cores);
    cout << endl;

    cout << "Done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs." << endl << endl;

    // Dumping pairwise result
    phmap::BinaryOutputArchive ar_out_1(string(input_prefix + "_pairwise." + part_name + ".bin").c_str());
    ar_out_1.saveBinary(edges->size());
    for (auto& [k, v] : *edges)
    {
        ar_out_1.saveBinary(k);
        ar_out_1.saveBinary(v);
    }

}