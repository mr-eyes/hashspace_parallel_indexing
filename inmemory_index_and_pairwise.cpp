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
#include <glob.h>
#include "parallel_hashmap/phmap_dump.h"
#include "lib.hpp"
#include "progressbar.hpp"

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;
typedef std::chrono::high_resolution_clock Time;


void pairwise(PAIRS_COUNTER* edges,
    flat_hash_map<uint64_t, std::vector<uint32_t>>* color_to_ids,
    flat_hash_map<uint64_t, uint32_t>& colorsCount,
    flat_hash_map<string, uint32_t>& groupName_to_kmerCount,
    flat_hash_map<uint64_t, string>& namesmap,
    int user_threads,
    string index_prefix);

void ranged_index(
    const vector<file_info>& all_files,
    flat_hash_map<uint64_t, std::vector<uint32_t>>* legend,
    flat_hash_map<uint64_t, uint32_t>& colorsCount,
    flat_hash_map<uint64_t, string>& namesmap,
    const int& kSize,
    const uint64_t& from_hash,
    const uint64_t& to_hash);

void inmemory_ranged_index(
    phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, phmap::flat_hash_set<uint64_t>>>,
    1,
    std::mutex>* loaded_bins,
    flat_hash_map<uint64_t, std::vector<uint32_t>>* legend,
    flat_hash_map<uint64_t, uint32_t>& colorsCount,
    flat_hash_map<uint64_t, string>& kSpider_namesmap,
    const int& kSize,
    const uint64_t& from_hash,
    const uint64_t& to_hash);


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
    int loading_cores = pairwise_cores;

    cout << "indexing cores: " << indexing_cores << endl;
    cout << "pairwise cores: " << pairwise_cores << endl;


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

    auto all_files = glob(bins_dir + "/*");
    bin_to_hashes->reserve(all_files.size());

    cout << "Loading all files and counting hashing using " << cores << " cores..." << endl;

    // Loading all bins
    progressbar bar(all_files.size());
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");
#pragma omp parallel num_threads(cores)
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

#pragma omp parallel num_threads(indexing_cores)
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

            auto* edges = new PAIRS_COUNTER();
            pairwise(
                edges,
                legend,
                colorsCount,
                groupName_to_kmerCount,
                namesmap,
                pairwise_cores,
                out_dir + "/part_" + to_string(i + 1)
            );

            delete legend;
            delete edges;

            /*
                     __   __        ___
                     |  \ /  \ |\ | |__
                     |__/ \__/ | \| |___

             */

             // TODO Merging pairwise



        }
    }
    delete bin_to_hashes;
}

