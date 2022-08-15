#include <iostream>
#include <cstdint>
#include <chrono>
#include <glob.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include "lib.hpp"

using namespace std;

uint64_t to_uint64_t(std::string const& value) {
    uint64_t result = 0;
    char const* p = value.c_str();
    char const* q = p + value.size();
    while (p < q) {
        result *= 10;
        result += *(p++) - '0';
    }
    return result;
}

vector<tuple<uint64_t, uint64_t>> splitted_ranges(uint64_t max_hash, int cores) {
    vector<tuple<uint64_t, uint64_t>> ranges;
    uint64_t from_hash = 0;
    uint64_t to_hash = 0;
    uint64_t step = (uint64_t)(max_hash / cores);
    for (int i = 0; i < cores; i++) {
        to_hash += step;
        ranges.push_back({ from_hash, to_hash });
        from_hash += step;
    }

    if (to_hash < max_hash)
        ranges[cores - 1] = { get<0>(ranges[cores - 1]), max_hash };

    return ranges;
}

vector<file_info> fetch(const std::string& pattern) {

    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        std::stringstream ss;
        ss << "glob() failed with return_value " << return_value << std::endl;
        throw std::runtime_error(ss.str());
    }

    std::vector<file_info> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {

        file_info finfo;

        finfo.path = std::string(glob_result.gl_pathv[i]);
        size_t _lastindex = finfo.path.find_last_of(".");
        finfo.prefix = finfo.path.substr(0, _lastindex);
        finfo.base_name = finfo.prefix.substr(finfo.prefix.find_last_of("/\\") + 1);

        filenames.push_back(finfo);
    }
    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}


void load_all_bins(string bins_dir, BINS_PHMAP* bin_to_hashes, BINS_KMER_COUNT* groupName_to_kmerCount, int cores) {
    auto all_files = fetch(bins_dir + "/*.bin");
    bin_to_hashes->reserve(all_files.size());
    groupName_to_kmerCount->reserve(all_files.size());
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

            groupName_to_kmerCount->try_emplace_l(_basename,
                [](BINS_KMER_COUNT::value_type& v) {},
                bin_size
            );

            bin_to_hashes->try_emplace_l(_basename,
                [](BINS_PHMAP::value_type& v) {},
                bin_hashes
            );

#pragma omp critical
            bar.update();
        }
    }
}
