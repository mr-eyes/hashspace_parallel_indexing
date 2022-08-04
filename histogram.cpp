#include <parallel_hashmap/phmap.h>
#include <vector>
#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include <glob.h>
#include <string>
#include <sstream>
#include <unordered_map>
#include "parallel_hashmap/phmap_dump.h"
#include<omp.h>

using namespace std;

using PAIRS_COUNTER = phmap::parallel_flat_hash_map<std::pair<uint64_t, uint64_t>,
    float, boost::hash<pair<uint64_t, uint64_t>>,
    std::equal_to<std::pair<uint64_t, uint64_t>>,
    std::allocator<std::pair<const std::pair<uint64_t, uint64_t>, uint64_t>>,
    4, std::mutex>;

vector<tuple<uint64_t, uint64_t>> splitted_ranges(uint64_t max_hash, int cores) {
    vector<tuple<uint64_t, uint64_t>> ranges;
    uint64_t from_hash = 0;
    uint64_t to_hash = 0;
    uint64_t step = (uint64_t)(max_hash / cores);
    cout << "step: " << step << endl;
    for (int i = 0; i < cores; i++) {
        to_hash += step;
        ranges.push_back({ from_hash, to_hash });
        from_hash += step;
    }

    if (to_hash < max_hash) {
        ranges[cores - 1] = { get<0>(ranges[cores - 1]), max_hash };
    }


    return ranges;
}

struct file_info {
    string path;
    string prefix;
    string extension;
    string basename;
};


vector<file_info> glob(const std::string& pattern) {

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
        finfo.basename = finfo.prefix.substr(finfo.prefix.find_last_of("/\\") + 1);

        filenames.push_back(finfo);
    }
    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

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

int main(int argc, char** argv) {
    if (argc < 4) {
        cout << "./histogram <bins_dir> <min_scale> <cores>" << endl;
        exit(1);
    }


    string bins_dir = argv[1];
    uint64_t scale = to_uint64_t(argv[2]);
    int user_threads = stoi(argv[3]);
    auto all_files = glob(bins_dir + "/*.bin");

    // if(cores > all_files.size()){
    //     cout << "lowering cores to " << all_files.size() << endl;
    //     cores = all_files.size();
    // }

    uint64_t max_hash = UINT64_MAX / (uint64_t)scale;
    // uint64_t max_hash = scale;
    cout << "max_hash: " << to_string(max_hash) << endl;

    auto ranges = splitted_ranges(max_hash, user_threads);
    PAIRS_COUNTER counter;


    for (int i = 0; i < ranges.size(); i++) {
        uint64_t from_hash = get<0>(ranges[i]);
        uint64_t to_hash = get<1>(ranges[i]);
        counter.insert(make_pair(make_pair(from_hash, to_hash), 0));
    }

    int processed = 0;
    int thread_num, num_threads, start, end, vec_i;
    int n = all_files.size();
    omp_set_num_threads(user_threads);

#pragma omp parallel private(vec_i,thread_num,num_threads,start,end)
    {
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        start = thread_num * n / num_threads;
        end = (thread_num + 1) * n / num_threads;
        for (vec_i = start; vec_i != end; ++vec_i) {

            auto file = all_files[vec_i];
            // cout << ++processed << "/" << all_files.size() << ":" << file.path << endl;
            phmap::flat_hash_set<uint64_t> bin_hashes;
            phmap::BinaryInputArchive ar_in(file.path.c_str());
            bin_hashes.phmap_load(ar_in);

            for (const auto& hash : bin_hashes) {
                for (int i = 0; i < ranges.size(); i++) {
                    uint64_t from_hash = get<0>(ranges[i]);
                    uint64_t to_hash = get<1>(ranges[i]);

                    if (!(from_hash <= hash && hash < to_hash)) {
                        continue;
                    }
                    auto _p = make_pair(from_hash, to_hash);
                    counter.try_emplace_l(_p, [](PAIRS_COUNTER::value_type& v) { v.second++; });
                }
            }

        }
    }

    cout << "dumping ..." << endl;
    ofstream phmapDump("stats.tsv");

    uint64_t total = 0;
    vector<uint64_t> all_counts;

    phmapDump << "from" << '\t' << "to" << '\t' << "count" << endl;

    for (const auto& v : ranges) {
        auto _count = (uint64_t)counter[make_pair(get<0>(v), get<1>(v))];
        total += _count;
        all_counts.push_back(_count);
        phmapDump << get<0>(v) << '\t' << get<1>(v) << '\t' << to_string(_count) << endl;
    }

    phmapDump.close();
    cout << "total kmers: " << total << endl;
    uint64_t mean = (uint64_t)total / user_threads;
    cout << "mean: " << to_string(mean) << endl;
#include <cmath>
    uint64_t _std = 0;
    for (int i = 0; i < all_counts.size(); ++i) {
        _std += pow(all_counts[i] - mean, 2);
    }
    _std = sqrt(_std / all_counts.size());
    cout << "standard deviation: " << to_string(_std) << endl;

}