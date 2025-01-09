//
// Created by nandgate on 9/30/24.
//

//#ifndef SGUFP_SOLVER_CUT_H
//#define SGUFP_SOLVER_CUT_H
#pragma once

#include <unordered_map>
#include <map>
#include <tuple>
#include <memory>
#include <utility>

using namespace std;

typedef int8_t shi; // short ints :)

enum CutType {
	OPTIMALITY,
	FEASIBILITY
};

inline static int cantor(int a, int b){
	return (a+b+1)* (a+b)/2 + b;
}

/**
 * hash function for tuple of three integers, to be used as key
 * for cut coefficients map.
 */
struct cut_tuple_hash{
	size_t operator()(const tuple<int,int,int>& tuple1) const{
		int a = get<0>(tuple1);
		int b = get<1>(tuple1);
		int c = get<2>(tuple1);
		return cantor(a, cantor(b,c));
	}
};

/**
 * equality of two tuples of three integers, to be used in comparing two keys in
 * cut coefficients map.
 */
struct cut_tuple_equal{
	bool operator()(const tuple<int,int,int>& tuple1, const tuple<int,int,int>& tuple2) const {
		return get<0>(tuple1) == get<0>(tuple2) &&
		        get<1>(tuple1) == get<1>(tuple2) &&
				get<2>(tuple1) == get<2>(tuple2);
	}
};

struct unordered_map_hash {
	/*
	 * Hash value of the unordered_map<tuple<int,int,int>>, double>
	 */
	size_t operator()(const unordered_map<tuple<int,int,int>, double>& coeff) {
		size_t seed = 0;
		for (const auto& pair: coeff) {
			// get hash of tuple
			auto h1 = cut_tuple_hash()(pair.first);
			auto h2 = hash<double>{}(pair.second);

			// TODO; use h1 and h2 in seed somehow.
			// naive hashing (INFO highly inefficient)
		}
		return seed;
	}
};

// typedef unordered_map<tuple<int,int,int>, double, cut_tuple_hash, cut_tuple_equal> CutCoefficients;
typedef map<tuple<int,int,int>,double> CutCoefficients; // LATER: change to unordered_map
class Cut{
	size_t hash;
public:
	bool operator==(const Cut& cut2) const { // ASAP fix this
		return (this->cutType == cut2.cutType) && (this->RHS == cut2.RHS) &&
				(this->cutCoeff == cut2.cutCoeff);
	}

	bool operator!=(const Cut& cut2) const {
		return *this == cut2;
	}

	size_t getHash(){return hash;}

	CutType cutType;
	double RHS;
	CutCoefficients cutCoeff;

	[[nodiscard]] double get(uint a, uint b, uint c) const {
		return cutCoeff.at(make_tuple(a,b,c));
	}

	Cut(CutType cutType_, double RHS_, CutCoefficients cutCoeff_):
		cutType{cutType_}, RHS{RHS_}, cutCoeff{std::move(cutCoeff_)}{
			// compute hash here.
	}
};

/**
 * hash function for Cut object.
 */
struct cut_hash{
	size_t operator()(const Cut& cut1) const {
		// bit manipulation of hash of RHS, cuttype, and coefficients.
		size_t h1 = hash<double>{}(cut1.RHS);
		size_t h2 = hash<int>{}(cut1.cutType);
		size_t h3 = h1 ^ (h2 << 2);
		// ASAP add hash function for cut coefficients.
		// DO NOT USE THIS FUNCTION WITHOUT THE HASH OF CUT COEFFICIENTS.
		size_t h4 = hash<int>{}(0);
		return h4 ^ (h3 << 1);
	}
};


static inline vector<vector<vector<shi>>> w2y(const vector<int>& w_solution, const shared_ptr<Network>& networkPtr){
	vector<shi> y_1 (networkPtr->n,0);
	vector<vector<shi>> y_2(networkPtr->n, y_1);
	vector<vector<vector<shi>>> y_bar(networkPtr->n, y_2);

	for (uint a = 0; a < w_solution.size(); a++){
		if (w_solution[a] != -1){
			auto arcId = networkPtr->processingOrder[a].second;
			auto q = networkPtr->networkArcs[arcId].headId;
			auto i = networkPtr->networkArcs[arcId].tailId;
			auto j = networkPtr->networkArcs[w_solution[a]].headId;
			y_bar[i][q][j] = 1;
			//cout << "i: " << i << ", q: " << q << ", j " << j << endl;
		}
	}
	return y_bar;
}

class CutContainer {
	// unordered_set of containers or vector of containers.
	#ifdef DEBUG
	void displayCutStats() const {
		string type = (cutType == FEASIBILITY) ? "FEASIBILITY" : "OPTIMALITY";
		cout << "********************** Cut stats for nerds ************************" << endl;
		cout << "Number of " << type <<" cuts: " << cuts.size() << endl;
		cout << "*******************************************************************" << endl;
	}
	#endif

public:
	vector<Cut> cuts;
	CutType cutType;

public:

	 // begin(){ return cuts.begin();};
	explicit CutContainer(CutType type_): cutType(type_){}

	bool empty() const noexcept {return cuts.empty();}
	bool isCutExists(const Cut& cut) {
	 	for (const auto& c : cuts) {
	 		if (c == cut) return true;
	 	}
	 	return false;
	}

	void insertCut(Cut cut) {
		cuts.push_back(cut);
	}

	void clearContainer() {
		cuts.clear();
	}

	~CutContainer() {

		#ifdef DEBUG
	 	// displayCutStats();
		#endif
		cuts.clear();
	}
};


namespace Inavap {

	struct q_word {
		uint16_t offset;
		uint16_t j;
		uint16_t i;
		uint16_t q;
	};

	static constexpr uint64_t Q_MASK = 0XFFFF;
	static constexpr uint64_t IQJ_MASK = 0XFFFFFFFFFFFF;
	static constexpr uint64_t IQ_MASK = 0XFFFFFFFF;

	class Cut {
		size_t hash_val;
		double RHS;
		// map<tuple<int,int,int>, double> coeff;
		vector<pair<uint64_t, double>> coeff;
		// vector<uint32_t> q_offsets; // The 16 MSB contains the offset and the 16 LSB is the 'q' element.

		// /**
		//  * Returns the offset corresponding to 'q' and 'i' in the cut. The function scans all the elements in the
		//  * offsets vector and returns the offset if match exists, else returns 0.
		//  */
		// [[nodiscard]] uint16_t getStart(uint16_t q) const noexcept {
		// 	// iterate through all the elements.
		// 	auto result = std::find_if(q_offsets.begin(), q_offsets.end(),
		// 		[&q](const auto p) { return !((p&IQ_MASK) ^ q);});
		// 	if (result != q_offsets.end()) return (*result)>>16; // return 16-MSB of the match.
		// 	return 0;
		// }

		/**
		 * Returns the start index corresponding to 'q' and 'i' in the cut. The function performs a jump-search on the
		 * i_q elements in the coeff vector if match exists, else returns 0.
		 */
		[[nodiscard]] uint16_t getStart(uint32_t qi) const noexcept {
			// jump-search: search all q_i with offset of current q_i.
			for (uint16_t i = 0; i < coeff.size(); ) {
				if ((coeff[i].first & IQ_MASK) ^ qi) i += coeff[i].first >> 48; // jump to next q_i.
				return i;
			}
			return 0;
		}

	public:

		explicit Cut(double RHS_, vector<pair<uint64_t, double>> coeff_): RHS{RHS_}, coeff{std::move(coeff_)} {
			hash_val = 0;
			/* hash function: sum of (index * key + hash(val)) */
			for (size_t i = 0; i< coeff.size(); i++) {
				size_t val_hash = std::hash<double>{}(coeff[i].second);
				size_t key_hash = coeff[i].first * i;
				hash_val += (key_hash ^ val_hash);
			}
		}

		explicit Cut(Cut&& c) noexcept :
				hash_val{move(c.hash_val)}, RHS{c.RHS} ,coeff{move(c.coeff)}{}

		bool operator==(const Cut& cut2) const {return cut2.hash_val == hash_val && cut2.RHS == RHS;}

		Cut& operator=(Cut&& c) noexcept {
			hash_val = c.hash_val;
			RHS = c.RHS;
			coeff = std::move(c.coeff);
			// q_offsets = move(c.q_offsets);
			return *this;
		}

		Cut& operator=(const Cut& cut2) {
			hash_val = cut2.hash_val;
			RHS = cut2.RHS;
			coeff = cut2.coeff;
			// q_offsets = cut2.q_offsets;
			return *this;
		}

		[[nodiscard]]size_t getHash() const noexcept {return hash_val;}

		/**
		 * Returns the cut coefficient for the given key if exists, returns 0 otherwise.
		 * The key is divided into 4 words (offset, j, i, q).
		 */
		[[nodiscard]] double get(uint64_t& key) const noexcept {

			/* The 16 MSB of the key might contain offset of the qi (from previous calls). If exists use it,
			 * else populate 16 MSB of key with the offset */
			uint64_t current = (key>>48);
			if (!current) {
				current = getStart((key & IQ_MASK)); // current now holds index (offset) of start of q_i in coeff.
				key |= (current<<48); // set offset in 16 MSB of key.
			}

			/* The reason for getStart() to return 0 if key didn't exist is, the conditional check in the below for loop
			 * will definitely fail in first iteration if the 'q' in the key doesn't match with 'q' in the element of cut. */

			for (; !((coeff[current].first & Q_MASK) ^ (key & Q_MASK)); ++current) {
				// extract iqj and compare with key.
				if (!((coeff[current].first & IQJ_MASK)^(key & IQJ_MASK))) return coeff[current].second;
			}
			return 0;
		}

		[[nodiscard]] double getRHS() const noexcept {return RHS;}

	};

	/**
	 * Returns the 64-bit value of bitwise 'or' of q,i,j at the respective bit positions.
	 */
	[[gnu::always_inline]] static uint64_t getKey(uint64_t q, uint64_t i, uint64_t j) {
		return (q | (i<<16) | (j<<32));
	}

	class CutContainer {
		vector<Cut> cuts;
	public:
		explicit CutContainer(size_t N = 128) {
			cuts.reserve(N);
		}

		CutContainer(CutContainer&& c) noexcept : cuts(std::move(c.cuts)){}

		CutContainer& operator=(CutContainer&& c) noexcept {
			cuts = std::move(c.cuts);
			return *this;
		}

		// TODO: define 'new' operator (efficient, without copying).
		void insertCut(Cut&& cut) {cuts.push_back(move(cut));}
		[[nodiscard]] bool isCutExists(const Cut& cut) const noexcept {return std::find(cuts.begin(), cuts.end(), cut) != cuts.end();}
		[[nodiscard]] size_t size() const noexcept {return cuts.size();}
		[[nodiscard]] bool empty() const noexcept {return cuts.empty();}
		void clearContainer() noexcept { cuts.clear();}

		auto begin() noexcept {return cuts.begin();}
		auto end() noexcept {return cuts.end();}
	};

	static Inavap::Cut cutToCut(const ::Cut& cut, const Network *networkPtr) {
		// Constructs the Inavap::Cut from ::Cut.
		vector<pair<uint64_t, double>> coeff;
		vector<uint32_t> q_offsets;
		// processing order defines the global order.
		for (const auto [index, arcId] : networkPtr->processingOrder) {
			const auto& arc = networkPtr->networkArcs[arcId];
			/* Not sure how the cast works on the bitwise operations.*/
			uint64_t i = arc.tailId;
			uint64_t q = arc.headId;
			uint offset = 0;
			for (uint64_t j: networkPtr->networkNodes[q].outNodeIds) {
				double val = cut.cutCoeff.at(make_tuple(i,q,j));
				if (val == 0.0) continue; // reduce size of cut.
				offset++;
				uint64_t key = Inavap::getKey(q,i,j);
				coeff.emplace_back(key,val);
			}
			if (!offset) continue; // set offset information in 16MSB of coeff.first
			uint64_t start = coeff.size() - offset;
			coeff[start].first |= (offset << 48);
		}
		return Cut{cut.RHS, coeff};
	}
}

//#endif //SGUFP_SOLVER_CUT_H
