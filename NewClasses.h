//
// Created by nandgate on 9/5/24.
//

#include <utility>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include "Network.h"
#include <iostream>
#include <algorithm>

using namespace std;

class Arc{
public:
	int id;
	int head;
	int tail;
	int value;
	int weight;

	Arc(){}

	Arc(int a, int b, int c, int d): id{a}, head{b}, tail{c}, value{d}, weight{0}{}
};

class Node{
public:
	int id;
	vector<int> parents;
	vector<int> children;
	unordered_set<int> states;
	int state2;

	Node(){};
	Node(int a): id{a}{}
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

class Diagram{

private:
	uint lastArc = 1;
	uint lastNode = 1;
	Type type;
	Prune strategy;
	uint maxWidth = 128;

	void reduceLayer(vector<uint>& currentLayer);
	void mergeNodes(Node& node1, Node& node2);

public:
	unordered_map<int,Node> nodes;
	unordered_map<int, Arc> arcs;
	vector<vector<int>> tree;
	int lastInserted = 0;
	Diagram() {

	}

	void build(const vector<pair<int,int>> processingOrder, Node& node, int index);
	//void build(const vector<pair<int,int>> processingOrder);

	vector<int>&& buildNextLayer(vector<int>& currentLayer, int index);
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
