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
#include <lib.hpp>
#include <ctime>
#include "parallel_hashmap/phmap_dump.h"

#define LOGGING 1

typedef std::chrono::high_resolution_clock Time;
using namespace std;


inline uint64_t to_uint64_t(std::string const& value) {
    uint64_t result = 0;
    char const* p = value.c_str();
    char const* q = p + value.size();
    while (p < q) {
        result *= 10;
        result += *(p++) - '0';
    }
    return result;
}


void bins_indexing(string bins_dir, string output_prefix, uint64_t from_hash, uint64_t to_hash) {

    kDataFrame* frame;
    std::string dir_prefix = bins_dir.substr(bins_dir.find_last_of("/\\") + 1);

    flat_hash_map<string, string> namesMap;
    string names_fileName = bins_dir;

    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;

    auto* legend = new phmap::parallel_flat_hash_map<uint64_t, std::vector<uint32_t>>();

    flat_hash_map<uint64_t, uint64_t> colorsCount;
    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;
    int detected_kSize = 0;

    int total_bins_number = 0;
    frame = new kDataFramePHMAP(31, mumur_hasher);

    flat_hash_map<string, string> basename_to_path;

    auto all_files = fetch(bins_dir + "/*");

    for (const auto& dirEntry : all_files) {
        string file_name = (string)dirEntry.path;
        size_t lastindex = file_name.find_last_of(".");
        string bin_prefix = file_name.substr(0, lastindex);
        std::string bin_basename = bin_prefix.substr(bin_prefix.find_last_of("/\\") + 1);


        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        if (extension != "bin") {
            cerr << "skipping " << file_name << " does not have extension .bin" << endl;
            continue;
        }

        basename_to_path.insert(pair(bin_basename, file_name));

        total_bins_number++;

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
            legend->insert(make_pair(groupID, tmp));
            colorsCount.insert(make_pair(groupID, 0));
            groupID++;
        }
    }


    // ----------------------------------------------------------------


    flat_hash_map<uint64_t, string> inv_groupNameMap;
    for (auto& _ : groupNameMap)
        inv_groupNameMap[_.second] = _.first;


    int currIndex = 0;
    string kmer;
    uint64_t tagBits = 0;
    uint64_t maxTagValue = (1ULL << tagBits) - 1;

    uint64_t lastTag = 0;
    readID = 0;

    int processed_bins_count = 0;
    auto begin_time = Time::now();
    uint_fast64_t current_kmers_numbers = 0;

    phmap::flat_hash_map<uint32_t, uint32_t> kmer_count;

    // START
    for (const auto& [bin_basename, bin_path] : basename_to_path) {
        //START
#if LOGGING == 1
        cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;
#endif
        flat_hash_map<uint64_t, uint64_t> convertMap;

        string readName = bin_basename;
        string groupName = bin_basename;

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

        begin_time = Time::now();
        phmap::flat_hash_set<uint64_t> bin_hashes;
        phmap::BinaryInputArchive ar_in(bin_path.c_str());
        bin_hashes.phmap_load(ar_in);

        // output kmer count
        kmer_count[groupNameMap[bin_basename]] = bin_hashes.size();


        for (const uint64_t& hashed_kmer : bin_hashes) {
            if (!(from_hash <= hashed_kmer && hashed_kmer < to_hash)) continue;
            uint64_t currentTag = frame->getCount(hashed_kmer);
            auto itc = convertMap.find(currentTag);
            if (itc == convertMap.end()) {
                vector<uint32_t> colors = legend->find(currentTag)->second;
                auto tmpiT = find(colors.begin(), colors.end(), readTag);
                if (tmpiT == colors.end()) {
                    colors.push_back(readTag);
                    sort(colors.begin(), colors.end());
                }

                string colorsString = to_string(colors[0]);
                for (int k = 1; k < colors.size(); k++) {
                    colorsString += ";" + to_string(colors[k]);
                }

                auto itTag = tagsMap.find(colorsString);
                if (itTag == tagsMap.end()) {
                    uint64_t newColor;
                    if (freeColors.size() == 0) {
                        newColor = groupID++;
                    }
                    else {
                        newColor = freeColors.top();
                        freeColors.pop();
                    }

                    tagsMap.insert(make_pair(colorsString, newColor));
                    legend->insert(make_pair(newColor, colors));
                    itTag = tagsMap.find(colorsString);
                    colorsCount[newColor] = 0;
                }
                uint64_t newColor = itTag->second;

                convertMap.insert(make_pair(currentTag, newColor));
                itc = convertMap.find(currentTag);
            }

            if (itc->second != currentTag) {

                colorsCount[currentTag]--;
                if (colorsCount[currentTag] == 0 && currentTag != 0) {

                    auto _invGroupNameIT = inv_groupNameMap.find(currentTag);
                    if (_invGroupNameIT == inv_groupNameMap.end()) {
                        freeColors.push(currentTag);
                        vector<uint32_t> colors = legend->find(currentTag)->second;
                        string colorsString = to_string(colors[0]);
                        for (int k = 1; k < colors.size(); k++) {
                            colorsString += ";" + to_string(colors[k]);
                        }
                        tagsMap.erase(colorsString);
                        legend->erase(currentTag);
                        if (convertMap.find(currentTag) != convertMap.end())
                            convertMap.erase(currentTag);
                    }

                }
                colorsCount[itc->second]++;
            }

            frame->setCount(hashed_kmer, itc->second);
        }
        readID += 1;
        groupCounter[groupName]--;
        if (colorsCount[readTag] == 0) {
            if (groupCounter[groupName] == 0) {
                freeColors.push(readTag);
                legend->erase(readTag);
            }

        }

#if LOGGING == 1
        auto loop_time_secs = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
        cout << "   loaded_kmers      " << bin_hashes.size() << endl;
        cout << "   uniq_added_kmers: " << frame->size() - current_kmers_numbers << endl;
        cout << "   total_kmers       " << frame->size() << " | load_factor: " << frame->load_factor() << endl;
        cout << "   total_colors      " << legend->size() << " | load_factor: " << legend->load_factor() << endl;
        cout << "   loop_time:        " << loop_time_secs << " secs" << endl;
        cout << "--------" << endl;
        current_kmers_numbers = frame->size();
#endif
        // END

    }
    delete frame;

#if LOGGING == 1
    cout << endl << endl;
    cout << "saving results ..." << endl;
    begin_time = Time::now();
    cout << "Indexing done!" << endl;
    cout << "dumping results ..." << endl;
#endif
    auto color_to_sources = new phmap::flat_hash_map<uint64_t, phmap::flat_hash_set<uint32_t>>();
    for (auto it : *legend) {
        phmap::flat_hash_set<uint32_t> tmp(std::make_move_iterator(it.second.begin()), std::make_move_iterator(it.second.end()));
        color_to_sources->operator[](it.first) = tmp;
    }

    phmap::BinaryOutputArchive ar_out_1(string(output_prefix + "_color_to_sources.bin").c_str());
    ar_out_1.saveBinary(color_to_sources->size());
    for (auto& [k, v] : *color_to_sources)
    {
        ar_out_1.saveBinary(k);
        ar_out_1.saveBinary(v);
    }

    phmap::BinaryOutputArchive ar_out_2(string(output_prefix + "_kmer_count.bin").c_str());
    kmer_count.phmap_dump(ar_out_2);

    phmap::BinaryOutputArchive ar_out_3(string(output_prefix + "_color_count.bin").c_str());
    colorsCount.phmap_dump(ar_out_3);

#if LOGGING == 1
    auto dumping_time = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
    cout << "Done saving results with prefix (" << output_prefix << ") in " << dumping_time << " secs." << endl;
#endif
}


int main(int argc, char** argv) {

    if (argc < 6) {
        cout << "run: ./ranged_index <bins_dir> <output_prefix> <this_part> <total_parts> <scale>" << endl;
        exit(1);
    }

    string bins_dir = argv[1];
    string output_prefix = argv[2];
    uint64_t this_part = to_uint64_t(argv[3]);
    uint64_t total_parts = to_uint64_t(argv[4]);
    uint64_t scale = to_uint64_t(argv[5]);


    if (this_part <= 0) {
        cout << "this_part must be >= 1" << endl;
        exit(1);
    }

    uint64_t max_hash = UINT64_MAX / (uint64_t)scale;
    auto ranges = splitted_ranges(max_hash, total_parts);
    uint64_t from_hash = get<0>(ranges[this_part - 1]);
    uint64_t to_hash = get<1>(ranges[this_part - 1]);

#if LOGGING == 1
    cout << "This run: " << endl;
    cout << "\t\t- from hash: " << to_string(from_hash) << endl;
    cout << "\t\t- to hash: " << to_string(to_hash) << endl;
    cout << "--------------------------" << endl;
#endif

    auto begin_time = Time::now();
    bins_indexing(bins_dir, output_prefix, from_hash, to_hash);
#if LOGGING == 1
    auto total = std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000;
    cout << "======================\n DONE PROCESSING PART " << this_part << "/" << total_parts << " in " << total << " secs." << endl;
#endif
}