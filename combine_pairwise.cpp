#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <vector>
#include <iostream>
#include <lib.hpp>
#include <fstream>
#include "progressbar.hpp"

#define CALC_ANI 1

#if CALC_ANI

#include <Python.h>

class toANI {
public:
    PyObject* moduleMainString, * moduleMain, * func;

    toANI() {
        Py_Initialize();
        PyRun_SimpleString("from sourmash.distance_utils import containment_to_distance\n");
        moduleMainString = PyUnicode_FromString("__main__");
        moduleMain = PyImport_Import(moduleMainString);
        PyRun_SimpleString("def to_ani(a,b,c,d): return containment_to_distance(a, b, c, n_unique_kmers=d).ani");
        func = PyObject_GetAttrString(moduleMain, "to_ani");
    }

    float calculate(float containment, long kSize, long scale, long n_unique_kmers) {
        PyObject* args = PyTuple_Pack(4,
            PyFloat_FromDouble(containment),
            PyLong_FromLong(kSize),
            PyLong_FromLong(scale),
            PyLong_FromLong(n_unique_kmers)
        );
        PyObject* result = PyObject_CallObject(func, args);
        return PyFloat_AsDouble(result);
    }
};




#endif

using namespace std;


void load_edges(const std::string& filename, PAIRS_COUNTER* final_map)
{
    auto* map = new PAIRS_COUNTER();
    phmap::BinaryInputArchive ar_in(filename.c_str());
    size_t size;
    ar_in.loadBinary(&size);
    map->reserve(size);
    while (size--)
    {
        std::pair<uint32_t, uint32_t> k;
        uint64_t v;
        ar_in.loadBinary(&k);
        ar_in.loadBinary(&v);
        map->insert_or_assign(std::move(k), std::move(v));
    }

    for (auto& [seq_pair, count] : *map) {
        final_map->try_emplace_l(seq_pair,
            [count](PAIRS_COUNTER::value_type& v) { v.second += count; },           // called only when key was already present
            count
        );
    }

    delete map;
}

inline flat_hash_map<string, int> parse_metadata(string input_prefix) {
    flat_hash_map<string, int> metadata_map;
    string line;
    ifstream fin(string(input_prefix + ".metadata").c_str());
    while (std::getline(fin, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, ',')) tokens.push_back(token);
        metadata_map[tokens[0]] = stoi(tokens[1]);
    }
    fin.close();
    return metadata_map;
}

int main(int argc, char** argv) {

    cout << "[DEBUG] argc: " << argc << endl;

#if CALC_ANI
    int expected_argc = 4;
#else
    int expected_argc = 3;
#endif



    if (argc != expected_argc) {
#if CALC_ANI
        cerr << "./combine_pairwise <input_prefix> <cores> <kSize>\n";
#else
        cerr << "./combine_pairwise <input_prefix> <cores>\n";
#endif
        exit(1);
    }

    string input_prefix = argv[1];
    int cores = stoi(argv[2]);
#if CALC_ANI
    int kSize = stoi(argv[3]);
#endif

    auto metadata_map = parse_metadata(input_prefix);

    cout << "detected scale: " << metadata_map["scale"] << endl;
    cout << "detected total_parts: " << metadata_map["parts"] << endl;

    vector<string> filenames;
    for (int part_id = 1; part_id <= metadata_map["parts"]; part_id++)
        filenames.push_back(input_prefix + "_pairwise.part" + to_string(part_id) + ".bin");

    if (cores > filenames.size()) {
        cout << "down scaling cores to " << filenames.size();
        cores = filenames.size();
    }

    flat_hash_map<uint32_t, uint32_t> group_id_to_kmer_count;
    {
        phmap::BinaryInputArchive ar_in(string(input_prefix + "_groupID_to_kmerCount.bin").c_str());
        group_id_to_kmer_count.phmap_load(ar_in);
    }
    assert(group_id_to_kmer_count.size());


    auto* edges = new PAIRS_COUNTER();



    int thread_num_1, num_threads_1, start_1, end_1, vec_i_1;
    int n_1 = filenames.size();
    omp_set_num_threads(cores);

    progressbar totalbar(n_1);
    totalbar.set_todo_char(" ");
    totalbar.set_done_char("█");
    totalbar.set_opening_bracket_char("[");
    totalbar.set_closing_bracket_char("]");

#pragma omp parallel private(vec_i_1,thread_num_1,num_threads_1,start_1,end_1)
    {
        thread_num_1 = omp_get_thread_num();
        num_threads_1 = omp_get_num_threads();
        start_1 = thread_num_1 * n_1 / num_threads_1;
        end_1 = (thread_num_1 + 1) * n_1 / num_threads_1;

        for (vec_i_1 = start_1; vec_i_1 != end_1; ++vec_i_1) {
            auto& filename = filenames[vec_i_1];
            load_edges(filename, edges);
            totalbar.update();
        }
    }
    cout << endl;

    std::ofstream myfile;
#if CALC_ANI
    myfile.open(input_prefix + "_kSpider_pairwise_with_ani.tsv");
#else
    myfile.open(input_prefix + "_kSpider_pairwise.tsv");
#endif
    myfile
        << "bin_1"
        << '\tbin_2'
        << '\tshared_kmers'
        << '\tmax_containment'
        << '\tavg_containment'
#if CALC_ANI
        << '\tavg_ani'
#endif
        << '\n';

    uint64_t line_count = 0;

#if CALC_ANI
    toANI ANI;
#endif


    for (const auto& edge : *edges) {
        // Skipping shared_kmers < 5
        if (edge.second < 5) continue;
        float cont_1_in_2 = (float)edge.second / group_id_to_kmer_count[edge.first.second];
        float cont_2_in_1 = (float)edge.second / group_id_to_kmer_count[edge.first.first];
        float avg_containment = (cont_1_in_2 + cont_2_in_1) / 2.0;
        float max_containment = max(cont_1_in_2, cont_2_in_1);
#if CALC_ANI
        float ani_1_in_2 = ANI.calculate(cont_1_in_2, kSize, metadata_map["scale"], metadata_map["scale"] * group_id_to_kmer_count[edge.first.second]);
        float ani_2_in_1 = ANI.calculate(cont_2_in_1, kSize, metadata_map["scale"], metadata_map["scale"] * group_id_to_kmer_count[edge.first.first]);
        float avg_ani = (ani_1_in_2 + ani_2_in_1) / 2;
#endif

        myfile
            << edge.first.first
            << '\t'
            << edge.first.second
            << '\t'
            << edge.second
            << '\t'
            << max_containment
            << '\t'
            << avg_containment
#if CALC_ANI
            << '\t' << avg_ani
#endif
            << '\n';
    }
    myfile.close();
    cout << "loaded " << edges->size() << " edges" << endl;
    std::ofstream metadata(input_prefix + ".metadata", std::ios_base::app | std::ios_base::out);
    metadata << "edges," << to_string(edges->size()) << "\n";
}