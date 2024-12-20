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
	#define RELAXED_STRATEGY 1
#endif
#ifndef EXACT_STRATEGY
	#define EXACT_STRATEGY 1
#endif

#ifndef MAX_WIDTH
	#define MAX_WIDTH 128
#endif

#ifndef NUMBERS_RESERVE
	#define NUMBERS_RESERVE 512
#endif

#ifndef DEBUG
	// #define DEBUG 1
	#define NDEBUG
	#define STATIC static
#else
	#define STATIC
#endif

#include <cassert>

#include "Network.h"
#include "Cut.h"
#include <cstdlib>
#include <set>
#include <algorithm>
#include <ctime>
#include <utility>
#include <optional>
#include <set>

#define assertm(exp, msg) assert(((void)msg, exp))

using namespace std;

typedef vector<ulint> vulint;

class Node_t {
public:
    vi states;
    vi solutionVector;
    double lb; // change to ints?
    double ub;
    uint globalLayer;

    Node_t() = default;

    Node_t(vi states_, vi solutionVector_, double lb_, double ub_, uint globalLayer_):
        states{std::move(states_)}, solutionVector{std::move(solutionVector_)},
        lb{lb_}, ub{ub_}, globalLayer{globalLayer_}{}
} ;

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
	uint nodeLayer = 0;
	uint globalLayer = 0;
	vector<ulint> incomingArcs;
	vector<ulint> outgoingArcs;
	set<int> states;
	double state2;
	vector<int> solutionVector;
	int objVal = INT32_MAX;

	DDNode():id{0}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{} {};
	explicit DDNode(ulint a): id{a}, incomingArcs{}, outgoingArcs{}, states{}, state2{0}, solutionVector{}{}
	DDNode(ulint id, uint layer, vi states_, vi solutionVector_): id{0}, incomingArcs{}, outgoingArcs{}, globalLayer{layer},
					states{states_.begin(), states_.end()}, solutionVector(std::move(solutionVector_)), state2{0}{}

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
			return n++;
		}

		void setNext(ulint x){
			numbers.push_back(x);
		}
	};

	const shared_ptr<Network> networkPtr;
	Number number;
	Type type;
	bool isExact = true;
	bool isInFeasible = false;
	bool isTreeDeleted = false; // true if the entire tree is deleted while refinement.
	vector<Node_t> cutSet;
	// info below two variables should be updated during tree compilation.
	uint startTree = 0; // the start position of the subtree in the global tree.
	int exactLayer = 0; // the position of exact layer with respect to root of subtree.
	unordered_set<ulint> deletedNodeIds; // deleted node ids during refinement on a single node.

	void updateTree();
	[[nodiscard]] vi computePathForExactNode(ulint nodeId) const;
	[[nodiscard]] vector<Node_t> generateExactCutSet() const;

	bool buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, bool stateChangesNext);
	ulint createChild(DDNode& parent, int decision);
	void buildNextLayer2(vector<ulint>& currentLayer, vector<ulint>& nextLayer);
	void buildNextLayer3(vector<ulint>& currentLayer, vector<ulint>& nextLayer);
	void buildNextLayer4(vector<ulint>& currentLayer, vector<ulint>& nextLayer);
	void buildNextLayer5(vector<ulint>& currentLayer, vector<ulint>& nextLayer);
	void buildNextLayer6(vector<ulint>& currentLayer, vector<ulint>& nextLayer);

public:

	unordered_map<ulint,DDNode> nodes;
	unordered_map<ulint, DDArc> arcs;
	vector<vector<ulint>> tree; // layer corresponds to vector of node ids.

	explicit DD(const shared_ptr<Network>& networkPtr_): networkPtr{networkPtr_}, type{RESTRICTED}{}
	explicit DD(const shared_ptr<Network>& networkPtr_, const Type type_): networkPtr{networkPtr_}, type{type_}{}

	optional<vector<Node_t>> build(DDNode &node);

	/// refinement helper functions ///

	void reduceLayer(vector<ulint> &currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void duplicateNode(ulint id);
	inline void updateState(const vector<ulint> &currentLayer, const set<int> &states);
	inline DDNode duplicate(const DDNode& node);

	/// refinement functions ///

	void applyFeasibilityCutRestricted(const Network& network, const Cut& cut);
	void applyFeasibilityCutRelaxed(const Network& network, const Cut& cut);
	void applyOptimalityCutRestricted(const Cut &cut);
	void applyOptimalityCutRelaxed(const Cut &cut);
	void applyOptimalityCut(const Network& network, const Cut& cut);
	void refineTree(const Network& network, Cut cut);
	void applyFeasibilityCut(const Network& network, const Cut& cut);

	double applyOptimalityCutRestrictedLatest(const Cut &cut);
	bool applyFeasibilityCutRestrictedLatest(const Cut &cut);
	double applyOptimalityCutHeuristic(const Cut &cut);
	bool applyFeasibilityCutHeuristic(const Cut &cut);

	/// node deletion functions ///

	void deleteArcById(ulint id);
	void deleteNodeById(ulint id);
	void removeNode(ulint id, bool isBatch=false);
	void batchRemoveNodes(const vulint& ids);
	void bottomUpDelete(ulint id);
	void topDownDelete(ulint id);
	void deleteNode(DDNode& node);
	void deleteArc(DDNode& parentNode, DDArc& arc, DDNode& childNode);

	/// getter functions ///
	int getExactLayer() const { return exactLayer;}
	int getGlobalPosition() const { return startTree; }
	[[nodiscard]] vector<Node_t> getExactCutSet();
	vi solution();
	bool isTreeExact() const {return isExact;}


	#ifdef DEBUG

	void displayArcLabels() const noexcept{
		cout << "\n ************************** Arcs **********************************" << endl;

		for (const auto& layer: tree) {
			for (auto id: layer) {
				const auto& node = nodes.at(id);
				for (auto outer : node.outgoingArcs) {
					const auto& arc = arcs.at(outer);
					cout << arc.decision <<" ";
				} cout << " : ";
			}
			cout << endl;
		}
	}
	void displayStats() const {
		cout << "\n*********************** DD Stats for nerds ************************" << endl;
		string ddtype;
		if (type == RESTRICTED) ddtype = "RESTRICTED";
		else if (type == RELAXED) ddtype = "RELAXED";
		else {
            #if EXACT_STRATEGY == 1
                ddtype = "EXACT (State-Reduction)";
            #else
                ddtype = "EXACT";
            #endif
		}
		cout << "Type: " << ddtype;
		if (startTree) {
			cout << " , Global order: " << startTree << endl;
		}
		else cout << " , Global Tree" << endl;
		cout << "Position of root in global tree: " << startTree << endl;
		cout << "Number of layers in tree: " << tree.size() - 2 << " + (root + terminal) = " << tree.size() << endl;
		cout << "Number of nodes: " << nodes.size() << endl;
		cout << "Number of arcs: " << arcs.size() << endl;
		cout << "Index of exact layer: " << exactLayer << " (contains " << tree[exactLayer].size() << " nodes)" << endl;
		cout << "Size of each layer : "; for (const auto& layer: tree) cout << layer.size() << " "; cout << endl;
		cout << "Size of each arc layer: "; for (const auto& layer: tree) { int count = 0;
			for (auto id: layer) count += nodes.at(id).outgoingArcs.size(); cout << count << " ";
		} cout << endl;
		cout << "*******************************************************************\n" << endl;
	}
	#endif
};

inline DDNode node2DDNode(const Node_t& node) {
	DDNode ddnode{0};
	ddnode.solutionVector = node.solutionVector;
	ddnode.states = set(node.states.begin(), node.states.end());
	ddnode.globalLayer = node.globalLayer;
	return ddnode;
}

/* custom hash functions for tuple and set */

struct set_hash{
	/*
	 * Hash of the entire set is Bitwise XOR of hash of individual elements.
	 */
	size_t operator()(const set<int>& s) const {
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
	size_t operator()(const tuple<set<int>,int,int>& t) const {
		size_t h1 = set_hash{}(get<0>(t));
		size_t h2 = hash<int>{}(get<1>(t));
		return h1 ^ (h2 << 2);
	}
};

struct tuple_equal {
	/*
	 * Two tuples are equal if first two elements of the tuple are equal.
	 */
	bool operator()(const tuple<set<int>,int,int>& t1,
	                const tuple<set<int>,int,int>& t2) const {
		auto& first = get<0>(t1);
		auto& second = get<0>(t2);
		return first == second && get<1>(t1) == get<1>(t2);
	}
};

/**
 * Returns a list of m unique random numbers from the interval [0,n)
 *
 * uses Fischer - Yates algorithm
 * @param n - range [0,n)
 * @param m
 * @return
 */
STATIC inline vector<uint> getShuffledList(const size_t n, const size_t m){
	assertm(n > m, "n must be greater than m");
	size_t temp_m = m;

	vui shuffle(n);
	for (size_t i = 0; i < n; i++) shuffle[i] = i;
	srand(time(nullptr));
	size_t current = n-1;
	// select m numbers from the shuffle.
	while (temp_m--){
		auto val = rand()%current;
		// shuffle a[val] and a[current]
		auto temp = shuffle[current];
		shuffle[current] = shuffle[val];
		shuffle[val] = temp;
		current--;
	}

	vui result(m);
	for (size_t i = 0; i < m; i++) result[i] = shuffle[n-m+i];
	// sort result
	std::sort(result.begin(), result.end());
	return result;
}
//
// class RestrictedDD {
// private:
// 	const shared_ptr<Network> networkPtr;
// 	ulint number = 1;
//
// 	bool isExact = false;
// 	bool isTreeDeleted = false;
// 	bool isTreeBuilt = false;
//
// 	const uint MAXWIDTH;
//
// 	uint startTree = 0;
// 	// uint exactLayer = 0;
//
// 	vector<Node_t> cutset;
//
// 	void updateStates(const vui& currentLayer, const unordered_set<int>& states);
//
// 	[[nodiscard]] vui buildNextLayer(const vui& currentLayer, bool hasStateChanged, bool& isExact, uint& nextSize);
//
// 	[[nodiscard]] vector<Node_t> generateExactCutSet(uint exactLayer) const;
//
// 	[[nodiscard]] vi computePathForNode(uint nodeId) const;
//
// 	/// deletion functions ///
//
// 	// void topDownDelete(uint id);
// 	// void bottomUpDelete(uint id);
// 	// void removeNode(uint id, bool isBatch);
// 	void batchRemoveNodes(const vui& nodeIds);
//
// 	void updateTree();
//
// 	// refinement functions
//
//
// public:
// 	vector<vector<uint>> tree;
// 	unordered_map<uint, DDNode> nodes;
// 	unordered_map<uint, DDArc> arcs;
//
// 	RestrictedDD(const shared_ptr<Network>& networkPtr_, const uint mw): networkPtr{networkPtr_}, MAXWIDTH{mw}{}
//
// 	[[nodiscard]] optional<vector<Node_t>> compile(DDNode root);
// 	/// refinement functions ///
// 	[[nodiscard]] bool applyFeasibilityCut(const Cut& cut) noexcept;
// 	[[nodiscard]] double applyOptimalityCut(const Cut& cut) noexcept;
// 	[[nodiscard]] vi solution() const noexcept;
//
// 	void displayStats() const noexcept {
//
// 	};
// };
//
// class RelaxedDD {
// 	const shared_ptr<Network> networkPtr;
//
// 	uint number = 1;
//
//
// };