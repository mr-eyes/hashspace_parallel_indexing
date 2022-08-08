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