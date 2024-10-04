//
// Created by nandgate on 6/1/24.
//

//#ifndef SGUFP_SOLVER_DD_H
//#define SGUFP_SOLVER_DD_H
#pragma once
// default prune strategy is removing/merging trailing nodes.
#ifndef PRUNE
	#define PRUNE TRAIL
#endif
#ifndef RESTRICTED_STRATEGY
	#define RESTRICTED_STRATEGY TRAIL
#endif
#ifndef RELAXED_STRATEGY
	#define RELAXED_STRATEGY MERGE
#endif

#ifndef MAX_WIDTH
	#define MAX_WIDTH (1<<10)
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
	// id, head and tail should be unsigned long ints.
	ulint id;
	ulint head; // id of the head node of the arc.
	ulint tail;
	int decision; // it should store the solution.
	double weight;

	DDArc(): id{0}, tail{0}, head{0}, decision{0}, weight{0}{}

	DDArc(ulint id_, ulint tail_, ulint head_, int decision_):
			id{id_}, tail{tail_}, head{head_}, decision{decision_}, weight{0}{}
};

class DDNode{
public:
	ulint id;
	vector<ulint> incomingArcs;
	vector<ulint> outgoingArcs;
	unordered_set<int> states;
	int state2;
	vector<int> solutionVector;
	int objVal = INT32_MAX; // set to max int value.

	DDNode():id{0}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{} {};
	explicit DDNode(ulint a): id{a}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{}{}
	//DDNode(const DDNode& node): id{node.id}, incomingArcs{node.incomingArcs}{} // TODO, copy constructor used in duplication node.

	// copy assignment.

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

enum Prune{
	TRAIL,
	RANDOM
};

class DD{

private:
	int lastArc = 1;
	int lastNode = 1;
	Type type = RESTRICTED;
	//int startTree = 0; // this variable denotes the position where the tree starts in the global order.
	vi cutset{};
	//int exactLayer = 0; // represents which layer is exact layer.
	// temporary

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

public:
	Number number;
	void reduceLayer(vector<ulint> &currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void deleteArcById(ulint id);
	void deleteNodeById(ulint id);
	void duplicateNode(ulint id);
	inline void updateState(const vector<ulint> &currentLayer, const unordered_set<int> &states);
	inline DDNode duplicate(const DDNode& node);
	vi solution();
	vector<DDNode> getExactCutset();
	vi computePathForExactNode(ulint nodeId);
//public:
	unordered_map<ulint,DDNode> nodes;
	unordered_map<ulint, DDArc> arcs;
	vector<vector<ulint>> tree;
	//ulint lastInserted = 1; // 0 is reserved for root node.
	int startTree = 0; // INFO these two should be updated during tree compilation.
	int exactLayer = 0;
	void applyOptimalityCut(const Network& network, const Cut& cut);
	void refineTree(const Network& network, Cut cut);
	void applyFeasibilityCut(const Network& network, const Cut& cut);
	// LATER add Network pointer to the DD class. remove Network parameter in all the functions.
	void removeNode(ulint id);
	void bottomUpDelete(ulint id);
	void topDownDelete(ulint id);

	unordered_set<ulint> deletedNodeIds;

	DD() {

	}

	void build(const Network& network, DDNode& node, int index);
	//void build(const vector<pair<int,int>> processingOrder);

	bool buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, int index);
};


struct set_hash{
	/*
	 * Hash decision of the entire set is Bitwise XOR of hash of individual elements.
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
	 * Hash decision of tuple considers only the first two elements of the tuple.
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
