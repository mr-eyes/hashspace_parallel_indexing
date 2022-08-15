#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <vector>
#include <iostream>
#include <lib.hpp>
#include <fstream>

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

    auto metadata_map = parse_metadata(input_prefix);

    cout << "detected scale: " << metadata_map["scale"] << endl;
    cout << "detected total_parts: " << metadata_map["parts"] << endl;

    vector<string> filenames;
    for (int part_id = 1; part_id <= metadata_map["parts"]; part_id++)
        filenames.push_back(input_prefix + "_pairwise.part" + to_string(part_id) + ".bin");

    // for (int i = 0; i < argc - 1; i++)
    //     filenames.push_back(argv[i + 1]);

    auto* edges = new PAIRS_COUNTER();


    for (auto& filename : filenames) {
        cout << "Loading " << filename << endl;
        load_edges(filename, edges);
    }

    std::ofstream myfile;
    myfile.open(input_prefix + "_kSpider_pairwise.tsv");
    myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';
    uint64_t line_count = 0;
    for (const auto& edge : *edges) {
        float max_containment = 0.0;
        // float max_containment = (float)edge.second / min(groupName_to_kmerCount[namesmap[edge.first.first]], groupName_to_kmerCount[namesmap[edge.first.second]]);
        myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << max_containment << '\n';
    }
    myfile.close();

    cout << "loaded " << edges->size() << " edges" << endl;

}