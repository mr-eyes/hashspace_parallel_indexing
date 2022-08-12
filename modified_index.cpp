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
#include "progressbar.hpp"

using namespace std;



using BINS_MAP = phmap::parallel_flat_hash_map<std::string, phmap::flat_hash_set<uint64_t>,
    phmap::priv::hash_default_hash<std::string>,
    phmap::priv::hash_default_eq<std::string>,
    std::allocator<std::pair<const std::string, phmap::flat_hash_set<uint64_t>>>,
    12,
    std::mutex
>;
using LEGENDS_MAP = phmap::parallel_flat_hash_map<uint64_t,
    std::vector<uint32_t>,
    std::hash<uint64_t>,
    std::equal_to<uint64_t>,
    std::allocator<std::pair<const uint64_t, vector<uint32_t>>>,
    4>; // 6 submaps because colors will grow

using LEGENDS_MAP_OLD = phmap::parallel_flat_hash_map<uint64_t, std::vector<uint32_t>>;





void bins_indexing_hashsplit(colored_kDataFrame* res, string bins_dir, int kSize, uint64_t from_hash, uint64_t to_hash) {

    kDataFrame* frame;
    std::string dir_prefix = bins_dir.substr(bins_dir.find_last_of("/\\") + 1);

    flat_hash_map<string, string> namesMap;
    string names_fileName = bins_dir;

    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;

    auto* legend = new LEGENDS_MAP();

    flat_hash_map<uint64_t, uint32_t> colorsCount;
    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;
    int detected_kSize = 0;

    int total_bins_number = 0;
    frame = new kDataFramePHMAP(kSize, mumur_hasher);

    flat_hash_map<string, string> basename_to_path;

    for (const auto& dirEntry : fetch(bins_dir + "/*")) {
        string file_name = dirEntry.path;
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

        basename_to_path.insert(make_pair(bin_basename, file_name));

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

    // // cout << "namesmap construction done..." << endl;


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



    // START
    for (const auto& [bin_basename, bin_path] : basename_to_path) {
        //START

        // cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

        flat_hash_map<uint64_t, uint64_t> convertMap;

        string readName = bin_basename;
        string groupName = bin_basename;

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

        phmap::flat_hash_set<uint64_t> bin_hashes;
        phmap::BinaryInputArchive ar_in(bin_path.c_str());
        bin_hashes.phmap_load(ar_in);

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
            // no need now
            // if (frame->getCount(hashed_kmer) != itc->second) {
            //     //frame->setC(kmer,itc->second);
            //     // cout << "Error Founded " << hashed_kmer << " from sequence " << readName << " expected "
            //         << itc->second << " found " << frame->getCount(hashed_kmer) << endl;
            //     exit(1);
            // }
        }
        readID += 1;
        groupCounter[groupName]--;
        if (colorsCount[readTag] == 0) {
            if (groupCounter[groupName] == 0) {
                freeColors.push(readTag);
                legend->erase(readTag);
            }

        }
        // cout << "   saved_kmers  (~" << frame->size() << ")." << endl << endl;
        // cout << "   saved_colors (~" << legend->size() << ")." << endl << endl;

        // END

    }


    colorTable* colors = new intVectorsTable();
    for (auto it : *legend) {
        colors->setColor(it.first, it.second);
    }

    // colored_kDataFrame* res = new colored_kDataFrame();
    res->setColorTable(colors);
    res->setkDataFrame(frame);
    for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
        uint32_t sampleID = groupNameMap[iit->second];
        res->namesMap[sampleID] = iit->second;
        res->namesMapInv[iit->second] = sampleID;
    }
    // cout << "saving to " << dir_prefix << " ..." << endl;
    // res->save(dir_prefix);
}



void ranged_index(
    const vector<file_info>& all_files,
    flat_hash_map<uint64_t, std::vector<uint32_t>>* legend,
    flat_hash_map<uint64_t, uint32_t>& colorsCount,
    flat_hash_map<uint64_t, string>& kSpider_namesmap,
    const int& kSize,
    const uint64_t& from_hash,
    const uint64_t& to_hash) {

    kDataFrame* frame;

    flat_hash_map<string, string> namesMap;
    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;



    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;
    int detected_kSize = 0;

    int total_bins_number = 0;
    frame = new kDataFramePHMAP(kSize, mumur_hasher);

    flat_hash_map<string, string> basename_to_path;

    for (const auto& file : all_files) {
        string file_name = file.path;
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

        basename_to_path.insert(make_pair(bin_basename, file_name));

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

    // // cout << "namesmap construction done..." << endl;


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



    // START
    for (const auto& [bin_basename, bin_path] : basename_to_path) {
        //START

        // cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

        flat_hash_map<uint64_t, uint64_t> convertMap;

        string readName = bin_basename;
        string groupName = bin_basename;

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

        phmap::flat_hash_set<uint64_t> bin_hashes;
        phmap::BinaryInputArchive ar_in(bin_path.c_str());
        bin_hashes.phmap_load(ar_in);

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

    }
    for (const auto& [grpID, grpName] : inv_groupNameMap) {
        kSpider_namesmap[grpID] = grpName;
    }
    delete frame;
}

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
    const uint64_t& to_hash) {

    kDataFrame* frame;

    flat_hash_map<string, string> namesMap;
    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;



    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;
    int detected_kSize = 0;

    int total_bins_number = 0;
    frame = new kDataFramePHMAP(kSize, mumur_hasher);

    flat_hash_map<string, string> basename_to_path;

    for (auto& __it : *loaded_bins) {

        string bin_basename = __it.first;

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

    // // cout << "namesmap construction done..." << endl;


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



    // START
    for (auto& [bin_basename, bin_hashes] : *loaded_bins) {
        //START

        // cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

        flat_hash_map<uint64_t, uint64_t> convertMap;

        string readName = bin_basename;
        string groupName = bin_basename;

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

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

    }
    for (const auto& [grpID, grpName] : inv_groupNameMap) {
        kSpider_namesmap[grpID] = grpName;
    }
    delete frame;
}




// ---------------------------------------------


void inmemory_ranged_index_single_for_kSPider(
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
    const uint64_t& to_hash,
    string output_prefix
) {

    kDataFrame* frame;

    flat_hash_map<string, string> namesMap;
    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;



    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;
    priority_queue<uint64_t, vector<uint64_t>, std::greater<uint64_t>> freeColors;
    flat_hash_map<string, uint64_t> groupCounter;
    int detected_kSize = 0;

    int total_bins_number = 0;
    frame = new kDataFramePHMAP(kSize, mumur_hasher);

    flat_hash_map<string, string> basename_to_path;

    for (auto& __it : *loaded_bins) {

        string bin_basename = __it.first;

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

    // // cout << "namesmap construction done..." << endl;


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

    // START
    // Loading all bins
    progressbar bar(inv_groupNameMap.size());
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");

    // START
    for (auto& [bin_basename, bin_hashes] : *loaded_bins) {
        //START

        // cout << "Processing " << ++processed_bins_count << "/" << total_bins_number << " | " << bin_basename << " ... " << endl;

        flat_hash_map<uint64_t, uint64_t> convertMap;

        string readName = bin_basename;
        string groupName = bin_basename;

        uint64_t readTag = groupNameMap.find(groupName)->second;


        convertMap.clear();
        convertMap.insert(make_pair(0, readTag));
        convertMap.insert(make_pair(readTag, readTag));

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
        bar.update();
    }
    cout << endl;

    for (const auto& [grpID, grpName] : inv_groupNameMap) {
        kSpider_namesmap[grpID] = grpName;
    }

    colorTable* colors = new intVectorsTable();
    for (auto it : *legend) {
        colors->setColor(it.first, it.second);
    }

    colored_kDataFrame* res = new colored_kDataFrame();
    res->setColorTable(colors);
    res->setkDataFrame(frame);
    for (auto iit = namesMap.begin(); iit != namesMap.end(); iit++) {
        uint32_t sampleID = groupNameMap[iit->second];
        res->namesMap[sampleID] = iit->second;
        res->namesMapInv[iit->second] = sampleID;
    }
    cout << "saving to " << output_prefix + "_idx_" << " ..." << endl;
    res->save(output_prefix + "_idx_");

    delete frame;
    delete res;
}