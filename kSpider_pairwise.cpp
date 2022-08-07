#include <iostream>
#include <cstdint>
#include <chrono>
#include "parallel_hashmap/phmap.h"
#include "kDataFrame.hpp"
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/functional/hash.hpp>
#include <ctime>
#include<omp.h>


using boost::adaptors::transformed;
using boost::algorithm::join;
using namespace std;
using namespace phmap;


typedef std::chrono::high_resolution_clock Time;


using PAIRS_COUNTER = phmap::parallel_flat_hash_map<std::pair<uint32_t, uint32_t>, std::uint64_t, boost::hash<pair<uint32_t, uint32_t>>, std::equal_to<std::pair<uint32_t, uint32_t>>, std::allocator<std::pair<const std::pair<uint32_t, uint32_t>, uint32_t>>, 12, std::mutex>;

template <typename T>
void ascending(T& dFirst, T& dSecond)
{
    if (dFirst > dSecond)
        std::swap(dFirst, dSecond);
}

class Combo {

public:
    Combo() = default;

    std::vector<std::pair<uint32_t, uint32_t>> combs;

    void combinations(int n) {
        this->combs.clear();
        this->comb(n, this->r, this->arr);
    }

private:
    int* arr = new int[2];
    int r = 2;

    void comb(int n, int r, int* arr) {
        for (int i = n; i >= r; i--) {
            // choose the first element
            arr[r - 1] = i;
            if (r > 1) { // if still needs to choose
                // recursive into smaller problem
                comb(i - 1, r - 1, arr);

            }
            else {
                this->combs.emplace_back(std::make_pair(arr[0] - 1, arr[1] - 1));
            }
        }
    }

};


void pairwise(
    PAIRS_COUNTER* edges,
    flat_hash_map<uint64_t, std::vector<uint32_t>>* color_to_ids,
    flat_hash_map<uint64_t, uint32_t>& colorsCount,
    flat_hash_map<string, uint32_t>& groupName_to_kmerCount,
    flat_hash_map<uint64_t, string>& namesmap,
    int user_threads,
    string index_prefix
) {

    auto begin_time = Time::now();
    clock_t begin_detailed_pairwise_comb, begin_detailed_pairwise_edges, begin_detailed_pairwise_edges_insertion;
    float detailed_pairwise_comb = 0.0;
    float detailed_pairwise_edges = 0.0;
    float detailed_pairwise_edges_insertion = 0.0;

    // convert map to vec for parallelization purposes.
    auto begin_col_to_ids = color_to_ids->begin();
    auto end_col_to_ids = color_to_ids->end();

    auto vec_color_to_ids = vector<pair<uint64_t, vector<uint32_t>>>(begin_col_to_ids, end_col_to_ids);

    int thread_num, num_threads, start, end, vec_i;
    int n = vec_color_to_ids.size();

    omp_set_num_threads(user_threads);
    begin_time = Time::now();

#pragma omp parallel private(vec_i,thread_num,num_threads,start,end)
    {
        thread_num = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        start = thread_num * n / num_threads;
        end = (thread_num + 1) * n / num_threads;

        for (vec_i = start; vec_i != end; ++vec_i) {
            auto& item = vec_color_to_ids[vec_i];
            Combo combo = Combo();
            combo.combinations(item.second.size());
            for (uint64_t i = 0; i < combo.combs.size(); i++) {
                // for (auto const& seq_pair : combo.combs) {
                auto const& seq_pair = combo.combs[i];
                uint32_t _seq1 = item.second[seq_pair.first];
                uint32_t _seq2 = item.second[seq_pair.second];
                ascending(_seq1, _seq2);

                auto _p = make_pair(_seq1, _seq2);
                uint32_t ccount = colorsCount[item.first];

                // edges->insert_or_assign(_p, ccount);// += ccount;

                edges->try_emplace_l(_p,
                    [ccount](PAIRS_COUNTER::value_type& v) { v.second += ccount; },           // called only when key was already present
                    ccount
                );

                // edges->try_emplace_l(_p,
                //     [ccount](PAIRS_COUNTER::value_type& v) { v.second+= ccount; }, // called only when key was already present
                //     [_p, ccount](const PAIRS_COUNTER::constructor& ctor) {
                //         ctor(_p, ccount); }
                //     );

            }
        }
    }

    // cout << "pairwise hashmap construction: " << std::chrono::duration<double, std::milli>(Time::now() - begin_time).count() / 1000 << " secs" << endl;
    // cout << "writing pairwise matrix to " << index_prefix << "_kSpider_pairwise.tsv" << endl;

    std::ofstream myfile;
    myfile.open(index_prefix + "_kSpider_pairwise.tsv");
    myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';
    uint64_t line_count = 0;
    for (const auto& edge : *edges) {
        float max_containment = (float)edge.second / min(groupName_to_kmerCount[namesmap[edge.first.first]], groupName_to_kmerCount[namesmap[edge.first.second]]);
        myfile << edge.first.first << '\t' << edge.first.second << '\t' << edge.second << '\t' << max_containment << '\n';
    }
    myfile.close();
}