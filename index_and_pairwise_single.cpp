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

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;
typedef std::chrono::high_resolution_clock Time;


int main(int argc, char** argv) {


    if (argc < 7) {
        cout << "./index <bins_dir> <output_prefix> <this_part> <total_parts> <scale> <cores>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string output_prefix = argv[2];
    uint64_t this_part = to_uint64_t(argv[3]);
    uint64_t total_parts = to_uint64_t(argv[4]);
    uint64_t scale = to_uint64_t(argv[5]);
    int cores = stoi(argv[6]);
    int loading_cores = cores / 2;
    if (loading_cores > 64) loading_cores = 64;
    else if (cores == 1) loading_cores = 1;

    if (this_part <= 0) {
        cout << "this_part must be >= 1" << endl;
        exit(1);
    }

    uint64_t max_hash = UINT64_MAX / (uint64_t)scale;
    auto ranges = splitted_ranges(max_hash, total_parts);
    uint64_t from_hash = get<0>(ranges[this_part - 1]);
    uint64_t to_hash = get<1>(ranges[this_part - 1]);

    int kSize = 21;

    // Constructing namesmap file for just a single time.
    flat_hash_map<string, uint32_t> groupName_to_kmerCount;
    flat_hash_map<uint32_t, uint32_t> groupID_to_kmerCount;
    flat_hash_map<uint64_t, string> namesmap;


    parallel_flat_hash_map<std::string, uint32_t,
        phmap::priv::hash_default_hash<std::string>,
        phmap::priv::hash_default_eq<std::string>,
        std::allocator<std::pair<const std::string, uint32_t>>,
        1,
        std::mutex> tmp_groupName_to_kmerCount;


    auto* bin_to_hashes = new parallel_flat_hash_map<std::string, flat_hash_set<uint64_t>,
        phmap::priv::hash_default_hash<std::string>,
        phmap::priv::hash_default_eq<std::string>,
        std::allocator<std::pair<const std::string, flat_hash_set<uint64_t>>>,
        1,
        std::mutex>();

    flat_hash_map<uint32_t, uint32_t> kmer_count;


    int total_bins_number = 0;
    uint32_t groupID = 0;
    int processed = 0;

    auto all_files = fetch(bins_dir + "/*");
    bin_to_hashes->reserve(all_files.size());

    cout << "Loading all files and counting hashes using " << cores << " cores..." << endl;

    // Loading all bins
    progressbar bar(all_files.size());
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");
#pragma omp parallel num_threads(loading_cores)
    {
#pragma omp for
        for (int x = 0; x < all_files.size(); x++) {
            const auto& file = all_files[x];
            size_t lastindex = file.path.find_last_of(".");

            std::string::size_type idx;
            idx = file.path.rfind('.');
            std::string extension = "";

            if (idx != std::string::npos) extension = file.path.substr(idx + 1);
            if (extension != "bin") {
                cerr << "skipping " << file.path << " does not have extension .bin" << endl;
                continue;
            }

            phmap::flat_hash_set<uint64_t> bin_hashes;
            phmap::BinaryInputArchive ar_in(file.path.c_str());
            bin_hashes.phmap_load(ar_in);

            auto bin_size = bin_hashes.size();
            auto _basename = file.base_name;

            tmp_groupName_to_kmerCount.try_emplace_l(_basename,
                [](parallel_flat_hash_map<string, uint32_t>::value_type& v) {},
                bin_size
            );

            bin_to_hashes->try_emplace_l(_basename,
                [](parallel_flat_hash_map<string, flat_hash_set<uint64_t>>::value_type& v) {},
                bin_hashes
            );

#pragma omp critical
            bar.update();

        }
    }

    for (const auto& [k, v] : tmp_groupName_to_kmerCount)
        groupName_to_kmerCount[k] = v;

    tmp_groupName_to_kmerCount.clear();

    bool ___conversion___ = false;

    std::ofstream fstream_kmerCount;
    fstream_kmerCount.open(output_prefix + "_kmer_count.tsv");
    fstream_kmerCount << "ID\tname\tseq\tkmers\n";
    uint64_t counter = 0;
    for (const auto& item : groupName_to_kmerCount) {
        fstream_kmerCount << item.first << '\t' << item.second << '\n';
    }
    fstream_kmerCount.close();

    cout << "\nDONE" << endl;

    cout << "Now, just wait..." << endl;

    auto* edges = new PAIRS_COUNTER();


    flat_hash_map<string, string> namesMap;
    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;
    flat_hash_map<uint64_t, string> inv_groupNameMap;
    string seqName, groupName;
    flat_hash_map<string, uint64_t> groupCounter;

    for (auto& __it : *bin_to_hashes) {

        string bin_basename = __it.first;

        seqName = bin_basename;
        groupName = bin_basename;

        namesMap.insert(make_pair(seqName, groupName));
        auto it = groupNameMap.find(groupName);
        groupCounter[groupName]++;
        if (it == groupNameMap.end()) {
            groupNameMap.insert(make_pair(groupName, groupID));
            tagsMap.insert(make_pair(to_string(groupID), groupID));
            vector<uint32_t> tmp;
            tmp.clear();
            tmp.push_back(groupID);
            groupID++;
        }
    }

    // Loading all bins
    cout << endl;

    /*
            _ _  _ ___  ____ _  _ _ _  _ ____
            | |\ | |  \ |___  \/  | |\ | | __
            | | \| |__/ |___ _/\_ | | \| |__]
    */

    flat_hash_map<uint64_t, uint32_t> colorsCount;
    auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();

    cout << "Indexing ..." << endl;
    auto begin_time = Time::now();
    inmemory_ranged_index_single_for_kSPider(
        bin_to_hashes,
        legend,
        colorsCount,
        inv_groupNameMap,
        kSize,
        from_hash,
        to_hash,
        output_prefix + "/part_" + to_string(this_part) + "_");


    for (auto& [_ID, _NAME] : namesmap)
        groupID_to_kmerCount[_ID] = groupName_to_kmerCount[_NAME];

    cout << "Indexing done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;


    /*
            __   __        ___
            |  \ /  \ |\ | |__
            |__/ \__/ | \| |___

    */

    /*

    ██   ██ ███████ ██████  ██ ██████  ███████ ██████
    ██  ██  ██      ██   ██ ██ ██   ██ ██      ██   ██
    █████   ███████ ██████  ██ ██   ██ █████   ██████
    ██  ██       ██ ██      ██ ██   ██ ██      ██   ██
    ██   ██ ███████ ██      ██ ██████  ███████ ██   ██

    */

    cout << "performing pairwise comparisons on " << legend->size() << " colors ..." << endl;
    begin_time = Time::now();
    pairwise(
        edges,
        legend,
        colorsCount,
        groupName_to_kmerCount,
        inv_groupNameMap,
        cores,
        output_prefix + "/part_" + to_string(this_part) + "_"
    );
    delete legend;
    cout << "Pairwise comparisons done in " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;


    /*
             __   __        ___
             |  \ /  \ |\ | |__
             |__/ \__/ | \| |___

     */

    delete edges;

    delete bin_to_hashes;
}

