//
// Created by nandgate on 3/25/2025.
//

#ifndef BENCH_CUT_H
#define BENCH_CUT_H

#include <iostream>
#include <benchmark/benchmark.h>
#include <algorithm>
#include <vector>
#include <random>
#include "../Cut.h"

using namespace std;


class NewCut {

	static constexpr uint64_t Q_MASK			= 0XFFFF;			// Turn on 16 LSB of the other operand.
	static constexpr uint64_t IQ_MASK			= 0XFFFFFFFF;		// Turn on 32 LSB of the other operand.
	static constexpr uint64_t IQJ_MASK			= 0XFFFFFFFFFFFF;	// Turn on 48 LSB of the other operand
	static constexpr uint64_t EXTRACT_16_LSB	= 0XFFFF;
	static constexpr uint64_t EXTRACT_I			= 0XFFFF0000;		// Turn on only the I bits of the other operand.
	static constexpr uint64_t EXTRACT_J			= 0XFFFF00000000;	// Turn on only the J bits of the other operand.

    vector<pair<uint64_t, double>> cutcoeffs;
    size_t hash;

    uint32_t getStart(uint32_t qi) {
        // return start value.
        for (uint16_t i = 0; i < cutcoeffs.size(); ) {
            // extract iq from both values and compare them.
            if ((cutcoeffs[i].first & IQ_MASK) ^ qi&IQ_MASK) {
                //auto offset = cutcoeffs[i].first>>48;
                //cout << "offset : " << offset << endl;
                i += (cutcoeffs[i].first>>48);
            }
            else { // cout << "Found match" << endl;
                return i;
            }
        }
        return 0;
    }

public:
    NewCut(const vector<pair<uint64_t, double>> &cutcoeffs, size_t hash) : cutcoeffs(cutcoeffs), hash(hash) {}

    double getLinear(uint64_t &key) {
        // do linear search and return value.
        for (const auto& [k,v] : cutcoeffs) {
            if (!(key&IQJ_MASK ^ k&IQJ_MASK)) return v;
        }
        return 0;
    }

    double get(uint64_t &key) {

        // if offset is empty, generate offset
        uint64_t current = key >> 48;
        // cout << "Existing offset: " << current << endl;
        // if (!current || true) {
            // fill current with offset)
            current = getStart(key&IQ_MASK);
            // key |= (current << 48);
        // }

        // iterate over all j's for valid i and q.
        for (; !((cutcoeffs[current].first & IQ_MASK) ^ (key & IQ_MASK)); ++current) {
            // extract iqj and compare with key.
            if (!((cutcoeffs[current].first & IQJ_MASK)^(key & IQJ_MASK))) return cutcoeffs[current].second;
        }
        return 0;
    }
};




static vector<uint64_t> getShuffledList(const size_t n, const size_t m){
	size_t temp_m = m;

	vector<uint64_t> shuffle(n);
	for (size_t i = 0; i < n; i++) shuffle[i] = i;
	srand(time(nullptr));
	size_t current = n-1;
	// select m numbers from the shuffle.
	while (temp_m-- && current){
		auto val = rand()%current;
		// shuffle a[val] and a[current]
		auto temp = shuffle[current];
		shuffle[current] = shuffle[val];
		shuffle[val] = temp;
		current--;
	}

	vector<uint64_t> result(m);
	for (size_t i = 0; i < m; i++) result[i] = shuffle[n-m+i];
	// sort result
	std::sort(result.begin(), result.end());
	return result;
}

static vector<pair<uint64_t, double>>  generateCut(size_t n) {
    // generate cuts

    std::random_device r;
	std::default_random_engine e(r());
	std::uniform_real_distribution<double> urd(-100,100);
	std::uniform_int_distribution<uint16_t> uid_q(0,25);
    std::uniform_int_distribution<uint16_t> uid_i(0,20);
    std::uniform_int_distribution<uint16_t> uid_j(0,10);

    vector<pair<uint64_t, double>> coeffs;
    for (int it = 0; it < n; ++it) {
        uint16_t q = uid_q(e);
        uint16_t i = uid_i(e);
        uint16_t j = uid_j(e);

        auto key = Inavap::getKey(q,i,j);
        double val = urd(e);
        coeffs.emplace_back(key, val);
    }

    // sort and fill offsets
    std::sort(coeffs.begin(), coeffs.end(), [](const auto& p1, const auto& p2) {
        const auto& k1 = p1.first;
        const auto& k2 = p2.first;
        // extract q and i from k1 and k2 and sort it
        auto qi1 = k1&Inavap::IQ_MASK;
        auto qi2 = k2&Inavap::IQ_MASK;
        return qi1 < qi2;
    });
    // fill offsets.
    unordered_map<uint32_t, pair<int, uint16_t>> counts;

    for (int i = 0; i < coeffs.size(); ++i) {
        const auto& key = coeffs[i].first;
        uint32_t k = key&Inavap::IQ_MASK;
        if (counts.find(k) == counts.end()) {
            counts[k] = {i,1};
        }
        else counts[k].second++;
    }
    for (const auto& [k,v] : counts) {
        auto [index,count] = v;
        uint64_t offset = static_cast<uint64_t>(count) << 48;
        coeffs[index].first |= offset;
    }

    return coeffs;
}

static vector<uint64_t> generateQueries(const vector<pair<uint64_t, double>>& coeffs, double percentage) {
    size_t n = coeffs.size();
    size_t m = n*percentage;
    auto indices = getShuffledList(n, m);

    std::sort(indices.begin(), indices.end());

    vector<uint64_t> queries;
    for (auto index : indices) {
        queries.push_back(coeffs[index].first);
    }
    return queries;
}
#endif //BENCH_CUT_H
