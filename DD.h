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
#ifndef MAX_WIDTH
	#define MAX_WIDTH 128
#endif


#include <utility>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <set>

#include "Network.h"
#include <iostream>
#include <algorithm>

using namespace std;

class DDArc{
public:
	int id;
	int head;
	int tail;
	int decision; // it should store the solution.
	int weight;

	DDArc(): id{0}, tail{0}, head{0}, decision{0}, weight{0}{}

	DDArc(int a, int b, int c, int d): id{a}, tail{b}, head{c}, decision{d}, weight{0}{}
};

class DDNode{
public:
	int id;
	vector<int> incomingArcs;
	vector<int> outgoingArcs;
	unordered_set<int> states;
	int state2;
	vector<int> solutionVector;
	int objVal = 0; // ASAP Update it during refinement.

	DDNode():id{0}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{} {};
	DDNode(int a): id{a}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{}{}
	//DDNode(const DDNode& node): id{node.id}, incomingArcs{node.incomingArcs}{} // ASAP, copy constructor used in duplication node.

	// copy assignment.

	~DDNode(){
		incomingArcs.clear();
		outgoingArcs.clear();
		states.clear();
		solutionVector.clear();
		// delete incomingArcs, outgoingArcs, states vectors.
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
	//Prune strategy;
	//int maxWidth = 128;
	//int startTree = 0; // this variable denotes the position where the tree starts in the global order.
	vi cutset{};
	//int exactLayer = 0; // represents which layer is exact layer.
	// temporary
	//vector<pair<int,int>> processingOrder;
	//unordered_map<int, vector<int>> stateUpdateMap;

public:
	void reduceLayer(vector<int> &currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void deleteArcById(int id);
	void deleteNodeById(int id);
	void duplicateNode(int id);
	inline void updateState(const vector<int> &currentLayer, const unordered_set<int> &states);
	inline DDNode duplicate(const DDNode& node);
	vi solution();
	vector<DDNode> getExactCutset();
	vi getSolutionVector(int nodeId);
//public:
	unordered_map<int,DDNode> nodes; // change to vector if needed.
	unordered_map<int, DDArc> arcs; // change to vector if needed.
	vector<vector<int>> tree;
	int lastInserted = 1; // 0 is reserved for root node.
	int startTree = 0;
	int exactLayer = 0;

	DD() {

	}

	void build(const Network& network, DDNode& node, int index);
	//void build(const vector<pair<int,int>> processingOrder);

	bool buildNextLayer(vector<int> &currentLayer, vector<int> &nextLayer, int index);
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
