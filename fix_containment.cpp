#include "parallel_hashmap/phmap.h"
#include "parallel_hashmap/phmap_dump.h"
#include <vector>
#include <iostream>
#include <lib.hpp>
#include <fstream>
#include "progressbar.hpp"

using namespace std;



int main(int argc, char** argv) {

    string input_prefix = argv[1];

    flat_hash_map<uint32_t, uint32_t> group_id_to_kmer_count;
    {
        phmap::BinaryInputArchive ar_in(string(input_prefix + "_groupID_to_kmerCount.bin").c_str());
        group_id_to_kmer_count.phmap_load(ar_in);
    }
    assert(group_id_to_kmer_count.size());



    std::ofstream myfile;
    myfile.open(input_prefix + "_kSpider_pairwise.tsv");
    myfile << "bin_1" << '\t' << "bin_2" << '\t' << "shared_kmers" << '\t' << "max_containment" << '\n';


    string line;
    ifstream fin(string(input_prefix + "_kSpider_pairwise.tsv.old").c_str());
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        std::vector<string> tokens;
        std::istringstream iss(line);
        std::string token;
        while (std::getline(iss, token, '\t')) tokens.push_back(token);
        uint32_t seq1 = stoi(tokens[0]);
        uint32_t seq2 = stoi(tokens[1]);
        uint32_t shared_kmers = stoi(tokens[2]);
        float max_containment = (float)shared_kmers / min(group_id_to_kmer_count[seq1], group_id_to_kmer_count[seq2]);
        myfile << seq1 << '\t' << seq2 << '\t' << shared_kmers << '\t' << max_containment << '\n';
    }

    fin.close();
}