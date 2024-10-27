//
// Created by nandgate on 9/30/24.
//

//#ifndef SGUFP_SOLVER_CUT_H
//#define SGUFP_SOLVER_CUT_H
#pragma once

#include <unordered_map>
#include <map>
#include <tuple>

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

// typedef unordered_map<tuple<int,int,int>, double, cut_tuple_hash, cut_tuple_equal> CutCoefficients;
typedef map<tuple<int,int,int>,double> CutCoefficients; // LATER: change to unordered_map
class Cut{
public:
	bool operator==(const Cut& cut2) const {
		return (this->cutType == cut2.cutType) && (this->RHS == cut2.RHS) &&
				(this->cutCoeff == cut2.cutCoeff);
	}

	bool operator!=(const Cut& cut2) const {
		return *this == cut2;
	}

	CutType cutType;
	double RHS;
	CutCoefficients cutCoeff;

	Cut(CutType cutType_, double RHS_, CutCoefficients cutCoeff_):
		cutType{cutType_}, RHS{RHS_}, cutCoeff{std::move(cutCoeff_)}{

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


static inline vector<vector<vector<shi>>> w2y(const vector<int>& w_solution, const Network& network){
	vector<shi> y_1 (network.n,0);
	vector<vector<shi>> y_2(network.n, y_1);
	vector<vector<vector<shi>>> y_bar(network.n, y_2);

	for (uint a = 0; a < w_solution.size(); a++){
		if (w_solution[a] != -1){
			auto arcId = network.processingOrder[a].second;
			auto q = network.networkArcs[arcId].headId;
			auto i = network.networkArcs[arcId].tailId;
			auto j = network.networkArcs[w_solution[a]].headId;
			y_bar[i][q][j] = 1;
			//cout << "i: " << i << ", q: " << q << ", j " << j << endl;
		}
	}
	return y_bar;
}

class CutContainer {
	// unordered_set of containers or vector of containers.
	vector<Cut> cuts;
	CutType cutType;

public:

	explicit CutContainer(CutType type_): cutType(type_){}

	bool isCutExists(const Cut& cut) {
		return true; // if cut exists? return true
	}

	void insertCut(Cut&& cut) {

	}

	void clearContainer() {
		cuts.clear();
	}

	~CutContainer() {
		cuts.clear();
	}
};

//#endif //SGUFP_SOLVER_CUT_H
