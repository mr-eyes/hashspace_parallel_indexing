#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <vector>
#include <iostream>
#include <lib.hpp>
#include <fstream>
#include "progressbar.hpp"

using namespace std;


void load_edges(const std::string& filename, PAIRS_COUNTER* final_map)
{
    auto* map = new PAIRS_COUNTER();
    phmap::BinaryInputArchive ar_in(filename.c_str());
    size_t size;
    ar_in.loadBinary(&size);
    map->reserve(size);
    while (size--)
    {
        std::pair<uint32_t, uint32_t> k;
        uint64_t v;
        ar_in.loadBinary(&k);
        ar_in.loadBinary(&v);
        map->insert_or_assign(std::move(k), std::move(v));
    }

    for (auto& [seq_pair, count] : *map) {
        final_map->try_emplace_l(seq_pair,
            [count](PAIRS_COUNTER::value_type& v) { v.second += count; },           // called only when key was already present
            count
        );
    }

    delete map;
}

inline flat_hash_map<string, int> parse_metadata(string input_prefix) {
    flat_hash_map<string, int> metadata_map;
    string line;
    ifstream fin(string(input_prefix + ".metadata").c_str());
    while (std::getline(fin, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, ',')) tokens.push_back(token);
        metadata_map[tokens[0]] = stoi(tokens[1]);
    }
    fin.close();
    return metadata_map;
}

int main(int argc, char** argv) {

    string input_prefix = argv[1];
    int cores = stoi(argv[2]);

    auto metadata_map = parse_metadata(input_prefix);

    cout << "detected scale: " << metadata_map["scale"] << endl;
    cout << "detected total_parts: " << metadata_map["parts"] << endl;

    vector<string> filenames;
    for (int part_id = 1; part_id <= metadata_map["parts"]; part_id++)
        filenames.push_back(input_prefix + "_pairwise.part" + to_string(part_id) + ".bin");

    if (cores > filenames.size()) {
        cout << "down scaling cores to " << filenames.size();
        cores = filenames.size();
    }

    flat_hash_map<uint32_t, uint32_t> group_id_to_kmer_count;
    {
        phmap::BinaryInputArchive ar_in(string(input_prefix + "_groupID_to_kmerCount.bin").c_str());
        group_id_to_kmer_count.phmap_load(ar_in);
    }
    assert(group_id_to_kmer_count.size());


    auto* edges = new PAIRS_COUNTER();



    int thread_num_1, num_threads_1, start_1, end_1, vec_i_1;
    int n_1 = filenames.size();
    omp_set_num_threads(cores);

    progressbar totalbar(n_1);
    totalbar.set_todo_char(" ");
    totalbar.set_done_char("█");
    totalbar.set_opening_bracket_char("[");
    totalbar.set_closing_bracket_char("]");

#pragma omp parallel private(vec_i_1,thread_num_1,num_threads_1,start_1,end_1)
    {
        thread_num_1 = omp_get_thread_num();
        num_threads_1 = omp_get_num_threads();
        start_1 = thread_num_1 * n_1 / num_threads_1;
        end_1 = (thread_num_1 + 1) * n_1 / num_threads_1;

        for (vec_i_1 = start_1; vec_i_1 != end_1; ++vec_i_1) {
            auto& filename = filenames[vec_i_1];
            load_edges(filename, edges);
            totalbar.update();
        }
    }
    cout << endl;

    std::ofstream myfile;
    myfile.open(input_prefix + "_kSpider_pairwise.tsv");
    myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\t' << "avg_containment" << '\n';
    uint64_t line_count = 0;
    for (const auto& edge : *edges) {
        float max_containment = (float)edge.second / min(group_id_to_kmer_count[edge.first.first], group_id_to_kmer_count[edge.first.second]);
        float min_containment = (float)edge.second / max(group_id_to_kmer_count[edge.first.first], group_id_to_kmer_count[edge.first.second]);
        float avg_containment = (max_containment + min_containment) / 2.0;

        myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << max_containment << '\t' << avg_containment << '\n';
    }
    myfile.close();

    cout << "loaded " << edges->size() << " edges" << endl;
    std::ofstream metadata(input_prefix + ".metadata", std::ios_base::app | std::ios_base::out);
    metadata << "edges," << to_string(edges->size()) << "\n";
}