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
	uint id;
	uint head;
	uint tail;
	int value;
	int weight;

	DDArc(){}

	DDArc(uint a, uint b, uint c, int d): id{a}, tail{b}, head{c}, value{d}, weight{0}{}
};

class DDNode{
public:
	int id;
	vector<int> incomingArcs;
	vector<int> outgoingArcs;
	unordered_set<int> states;
	int state2;
	vector<int> solutionVector;

	DDNode(){};
	DDNode(int a): id{a}{}
	DDNode(const DDNode& node){} // ASAP, constructor used in duplication node.

	// copy assignment.

	~DDNode(){
		// destructor.
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
	uint lastArc = 1;
	uint lastNode = 1;
	Type type;
	Prune strategy;
	uint maxWidth = 128;

	// temporary
	//vector<pair<int,int>> processingOrder;
	//unordered_map<int, vector<int>> stateUpdateMap;

	void reduceLayer(vector<uint>& currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void deleteArcById(uint id);
	void deleteNodeById(uint id);
	void duplicateNode(uint id);

public:
	unordered_map<uint,DDNode> nodes; // change to vector if needed.
	unordered_map<uint, DDArc> arcs; // change to vector if needed.
	vector<vector<uint>> tree;
	int lastInserted = 1; // 0 is reserved for root node.

	DD() {

	}

	void build(const Network& network, DDNode& node, int index);
	//void build(const vector<pair<int,int>> processingOrder);

	void buildNextLayer(vector<uint>& currentLayer, vector<uint>& nextLayer, int index);
};


struct set_hash{
	/*
	 * Hash value of the entire set is Bitwise XOR of hash of individual elements.
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
	 * Hash value of tuple considers only the first two elements of the tuple.
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
