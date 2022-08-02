#include "algorithms.hpp"
#include <parallel_hashmap/phmap.h>

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
    ranges[cores - 1] = { from_hash, max_hash };
    return ranges;
}

int main() {

    auto ranges = splitted_ranges(18446744073709551, 64);
    cout << "ranges size: ";
    cout << ranges.size() << endl;
    for(int i = 0; i < ranges.size(); i++){
        cout << i+1 <<"- (" << get<0>(ranges[i]) << ", " << get<1>(ranges[i]) << ")" << endl;
    }

}

