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

    if (argc < 4) {
        cout << "./index <bins_dir> <out_dir> <min_scale> <cores>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string out_dir = argv[2];
    int min_scale = stoi(argv[3]);
    int cores = stoi(argv[4]);
    int indexing_cores = (int)(0.9 * cores);
    int pairwise_cores = cores - indexing_cores;
    int loading_cores = (int)(cores * 0.5);

    // cout << "indexing cores: " << indexing_cores << endl;
    // cout << "pairwise cores: " << pairwise_cores << endl;


    uint64_t max_hash = UINT64_MAX / (uint64_t)min_scale;

    auto ranges = splitted_ranges(max_hash, cores);

    int kSize = 21;

    // Constructing namesmap file for just a single time.
    flat_hash_map<string, uint32_t> groupName_to_kmerCount;

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


    int total_bins_number = 0;
    uint32_t groupID = 0;
    int processed = 0;

    auto all_files = fetch(bins_dir + "/*");
    bin_to_hashes->reserve(all_files.size());

    cout << "Loading all files and counting hashing using " << cores << " cores..." << endl;

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




    std::ofstream fstream_kmerCount;
    fstream_kmerCount.open(out_dir + "/kSpider_kmer_count.tsv");
    fstream_kmerCount << "ID\tname\tseq\tkmers\n";
    uint64_t counter = 0;
    for (const auto& item : groupName_to_kmerCount) {
        fstream_kmerCount << item.first << '\t' << item.second << '\n';
    }
    fstream_kmerCount.close();

    cout << "DONE" << endl;

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

    for (auto& _ : groupNameMap)
        inv_groupNameMap[_.second] = _.first;


#pragma omp parallel num_threads(cores)
    {
#pragma omp for

        for (int i = 0; i < ranges.size(); i++) {
            uint64_t from_hash = get<0>(ranges[i]);
            uint64_t to_hash = get<1>(ranges[i]);

            // cout << "from hash: " << to_string(from_hash) << endl;
            // cout << "to hash: " << to_string(to_hash) << endl;
            // cout << "max hash " << to_string(max_hash) << endl;

            // cout << "INDEXING PART " << i + 1 << endl;

            /*
                    _ _  _ ___  ____ _  _ _ _  _ ____
                    | |\ | |  \ |___  \/  | |\ | | __
                    | | \| |__/ |___ _/\_ | | \| |__]
            */

            flat_hash_map<uint64_t, uint32_t> colorsCount;
            auto* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
            flat_hash_map<uint64_t, string> namesmap;







            inmemory_ranged_index(
                bin_to_hashes,
                legend,
                colorsCount,
                namesmap,
                kSize,
                from_hash,
                to_hash);









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

            pairwise_aggregate(
                edges,
                legend,
                colorsCount,
                groupName_to_kmerCount,
                inv_groupNameMap,
                1 // cores
            );

            delete legend;

            /*
                     __   __        ___
                     |  \ /  \ |\ | |__
                     |__/ \__/ | \| |___

             */

             // TODO Merging pairwise



        }
    }

    std::ofstream myfile;
    myfile.open(out_dir + "/kSpider_pairwise.tsv");
    myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';
    uint64_t line_count = 0;
    for (const auto& edge : *edges) {
        float max_containment = (float)edge.second / min(groupName_to_kmerCount[inv_groupNameMap[edge.first.first]], groupName_to_kmerCount[inv_groupNameMap[edge.first.second]]);
        // myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << max_containment << '\n';
        myfile << inv_groupNameMap[edge.first.first] << '\t' << inv_groupNameMap[edge.first.second] << '\t' << edge.second << '\t' << max_containment << '\n';

    }
    myfile.close();

    delete edges;

    delete bin_to_hashes;
}

