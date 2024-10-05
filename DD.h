//
// Created by nandgate on 6/1/24.
//

//#ifndef SGUFP_SOLVER_DD_H
//#define SGUFP_SOLVER_DD_H
#pragma once
// default prune strategy is removing/merging trailing nodes.
/*
 * RESTRICTED_STRATEGY
 *  1 - TRAIL
 *  2 - RANDOM
 *
 * RELAXED_STRATEGY
 *  1 - TRAIL
 *  2 - PARENT_CHILD
 */

#ifndef RESTRICTED_STRATEGY
	#define RESTRICTED_STRATEGY 1
#endif
#ifndef RELAXED_STRATEGY
	#define RELAXED_STRATEGY 2
#endif

#ifndef MAX_WIDTH
	#define MAX_WIDTH 128
#endif

#ifndef NUMBERS_RESERVE
	#define NUMBERS_RESERVE 512
#endif

#include "Network.h"
#include "Cut.h"

#include <set>
#include <algorithm>

using namespace std;

/*
 * INFO Ids of DDArc and DDNode are unsigned long ints.
 */

class DDArc{
public:
	ulint id; // key for arcs map.
	ulint head; // id of the outgoing node.
	ulint tail; // id of the incoming node.
	int decision; // decision of the variable
	double weight; // weight

	DDArc(): id{0}, head{0}, tail{0}, decision{0}, weight{0}{}

	DDArc(ulint id_, ulint tail_, ulint head_, int decision_):
			id{id_}, head{head_}, tail{tail_}, decision{decision_}, weight{0}{}
};

class DDNode{
public:
	ulint id;
	vector<ulint> incomingArcs;
	vector<ulint> outgoingArcs;
	unordered_set<int> states;
	double state2;
	vector<int> solutionVector;
	int objVal = INT32_MAX;

	DDNode():id{0}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{} {};
	explicit DDNode(ulint a): id{a}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{}{}

	~DDNode(){
		incomingArcs.clear();
		outgoingArcs.clear();
		states.clear();
		solutionVector.clear();
	}
};

enum Type {
	RESTRICTED,
	RELAXED,
	EXACT
};

class DD{
private:

	class Number{
	private:
		ulint n; // INFO n starts from 1. 0 is reserved for (sub-)root.
		vector<ulint> numbers;
	public:
		Number(): n{1}, numbers{}{
			numbers.reserve(NUMBERS_RESERVE);
		}

		ulint getNext(){
			if (!numbers.empty()) {
				uint x = numbers.back();
				numbers.pop_back();
				return x;
			}
			else return n++;
		}

		void setNext(ulint x){
			numbers.push_back(x);
		}
	};
	Type type;
	vi cutset{};

public:
	Number number;
	unordered_map<ulint,DDNode> nodes;
	unordered_map<ulint, DDArc> arcs;
	vector<vector<ulint>> tree; // layer corresponds to vector of node ids.
	// info below two variables should be updated during tree compilation.
	int startTree = 0; // the start position of the subtree in the global tree.
	int exactLayer = 0; // the position of exact layer with respect to root of subtree.
	unordered_set<ulint> deletedNodeIds; // deleted node ids during refinement on a single node.


	bool buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, int index);

	/// refinement helper functions ///

	void reduceLayer(vector<ulint> &currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void duplicateNode(ulint id);
	inline void updateState(const vector<ulint> &currentLayer, const unordered_set<int> &states);
	inline DDNode duplicate(const DDNode& node);

	/// refinement functions ///

	void applyFeasibilityCutRestricted(const Network& network, const Cut& cut);
	void applyOptimalityCut(const Network& network, const Cut& cut);
	void refineTree(const Network& network, Cut cut);
	void applyFeasibilityCut(const Network& network, const Cut& cut);
	// LATER add Network pointer to the DD class. remove Network parameter in all the functions.

	/// node deletion functions ///

	void deleteArcById(ulint id);
	void deleteNodeById(ulint id);
	void removeNode(ulint id);
	void bottomUpDelete(ulint id);
	void topDownDelete(ulint id);

	/// other functions ///

	vi computePathForExactNode(ulint nodeId);


	DD(): type{RESTRICTED}{}
	explicit DD(Type type_): type{type_}{}

	void build(const Network& network, DDNode& node, int index);
	vi solution();
	vector<DDNode> getExactCutset();

};

/* custom hash functions for tuple and set */

struct set_hash{
	/*
	 * Hash of the entire set is Bitwise XOR of hash of individual elements.
	 */
	size_t operator()(const unordered_set<int>& s) const {
		size_t h = 0;
		for (auto i: s)
			h ^= hash<int>{}(i);
		return h;
	}
};

struct tuple_hash{
	/*
	 * Hash of tuple is computed with only the first two elements of the tuple.
	 */
	size_t operator()(const tuple<unordered_set<int>,int,int>& t) const {
		size_t h1 = set_hash{}(get<0>(t));
		size_t h2 = hash<int>{}(get<1>(t));
		return h1 ^ (h2 << 2);
	}
};

struct tuple_equal {
	/*
	 * Two tuples are equal if first two elements of the tuple are equal.
	 */
	bool operator()(const tuple<unordered_set<int>,int,int>& t1,
	                const tuple<unordered_set<int>,int,int>& t2) const {
		auto& first = get<0>(t1);
		auto& second = get<0>(t2);
		return first == second && get<1>(t1) == get<1>(t2);
	}
};
