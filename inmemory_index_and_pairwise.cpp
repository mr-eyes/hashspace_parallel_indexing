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
// #include "argparse.hpp"

using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;
typedef std::chrono::high_resolution_clock Time;


int main(int argc, char** argv) {

    // argparse::ArgumentParser program("kSpider_variant");
    // program.add_argument("-b", "--bins-dir")
    //     .required()
    //     .help("bins directory (should exist).");
    // program.add_argument("-o", "--output-dir")
    //     .required()
    //     .help("output directory (should exist).");
    // program.add_argument("-s", "--scale")
    //     .required()
    //     .help("the minimum scale found in the bins")
    //     .scan<'i', int>();
    // program.add_argument("-c", "--cores")
    //     .required()
    //     .help("number of cores")
    //     .scan<'i', int>();
    // program.add_argument("-p", "--parts")
    //     .required()
    //     .help("parts index will split into")
    //     .scan<'i', int>();


    // try {
    //     program.parse_args(argc, argv);
    // }
    // catch (const std::runtime_error& err) {
    //     std::cerr << err.what() << std::endl;
    //     std::cerr << program;
    //     std::exit(1);
    // }


    // string bins_dir = program.get<std::string>("-b");
    // string out_dir = program.get<std::string>("-o");
    // int min_scale = program.get<int>("-s");
    // int cores = program.get<int>("-c");
    // int split_parts = program.get<int>("-p");
    // int indexing_cores = (int)(0.9 * cores);
    // int pairwise_cores = cores - indexing_cores;
    // int loading_cores = (int)(cores * 0.5);

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
    int split_parts = cores;



    uint64_t max_hash = UINT64_MAX / (uint64_t)min_scale;

    auto ranges = splitted_ranges(max_hash, split_parts);

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
    fstream_kmerCount.open(out_dir + "/kSpider_kmer_count.tsv");
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
    progressbar totalbar(cores * 2);
    totalbar.set_todo_char(" ");
    totalbar.set_done_char("█");
    totalbar.set_opening_bracket_char("[");
    totalbar.set_closing_bracket_char("]");

    int thread_num_1, num_threads_1, start_1, end_1, vec_i_1;
    int n_1 = ranges.size();
    omp_set_num_threads(cores);

#pragma omp parallel private(vec_i_1,thread_num_1,num_threads_1,start_1,end_1)
    {

        thread_num_1 = omp_get_thread_num();
        num_threads_1 = omp_get_num_threads();
        start_1 = thread_num_1 * n_1 / num_threads_1;
        end_1 = (thread_num_1 + 1) * n_1 / num_threads_1;

        for(vec_i_1 = start_1; vec_i_1 != end_1; ++vec_i_1) {
            uint64_t from_hash = get<0>(ranges[vec_i_1]);
            uint64_t to_hash = get<1>(ranges[vec_i_1]);

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






            inmemory_ranged_index(
                bin_to_hashes,
                legend,
                colorsCount,
                namesmap,
                kSize,
                from_hash,
                to_hash);




#pragma omp critical
            {
                totalbar.update();
                if (!___conversion___) {
                    for (auto& [_ID, _NAME] : namesmap) {
                        groupID_to_kmerCount[_ID] = groupName_to_kmerCount[_NAME];
                        ___conversion___ = true;
                    }
                }
            }

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


#pragma omp critical
            totalbar.update();

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
    for (auto& edge : *edges) {
        float max_containment = (float)edge.second / min(groupID_to_kmerCount[edge.first.first], groupID_to_kmerCount[edge.first.second]);
        myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << to_string(max_containment) << '\n';
        // myfile << namesmap[edge.first.first] << '\t' << namesmap[edge.first.second] << '\t' << edge.second << '\t' << max_containment << '\n';

    }
    myfile.close();

    delete edges;

    delete bin_to_hashes;
}

