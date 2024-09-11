//
// Created by nandgate on 6/1/24.
//

//#ifndef SGUFP_SOLVER_DD_H
//#define SGUFP_SOLVER_DD_H
#pragma once


#include <utility>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include "Network.h"
#include <iostream>
#include <algorithm>

using namespace std;

class DDArc{
public:
	int id;
	int head;
	int tail;
	int value;
	int weight;

	DDArc(){}

	DDArc(int a, int b, int c, int d): id{a}, head{b}, tail{c}, value{d}, weight{0}{}
};

class DDNode{
public:
	int id;
	vector<int> incomingArcs;
	vector<int> outgoingArcs;
	unordered_set<int> states;
	int state2;

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

	void reduceLayer(vector<uint>& currentLayer);
	void mergeNodes(DDNode& node1, DDNode& node2);
	void deleteArcById(uint id);
	void deleteNodeById(uint id);
	void duplicateNode(uint id);

public:
	unordered_map<uint,DDNode> nodes;
	unordered_map<uint, DDArc> arcs;
	vector<vector<uint>> tree;
	int lastInserted = 0;
	DD() {

	}

	void build(const vector<pair<int,int>> processingOrder, DDNode& node, int index);
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

//
//#include <iostream>
//#include <vector>
//#include <forward_list>
//#include <unordered_set>
//#include "Network.h"
//#include <execution>
//
//using namespace std;
//
//typedef vector<uint32_t> vui; // alias for vector of 32-bit unsigned integers.
//typedef vector<int> vi; // alias for vector of 32-bit signed integers.
//
//enum PruneStrategy {RANDOM /* removes the nodes at random in the given layer.*/ ,
//		TRAIL /*removes the trailing nodes that are past the max width. */};
//
//enum DDType {
//	RESTRICTED /* RESTRICTED DD */,
//	RELAXED,
//	EXACT
//};
//
//class DDArc{
//public:
//	uint id;
//	uint headId; // id of child node in the next layer of tree.
//	uint tailId;
//	int label; // decision.
//
//	DDArc(uint arc_id, uint head_id, uint tail_id, int arc_label): id{arc_id}, headId{head_id}, tailId{tail_id}{}
//
//	DDArc(uint arc_id, uint head_id, int arc_label)
//		: id{arc_id}, headId{head_id}, label{arc_label} {}
//};
//
//class DDNode {
//public:
//	uint nodeId; /* id of the node in the current layer. (assigned sequentially)  */
//	uint64_t objectiveVal; /* objective value is an unsigned long number */
//	forward_list<uint32_t> state; /* node might have multiple state values. */
//	uint layerNo; /* with respect to global ordering of the variables */
//	bool isExact;
//	vector<uint> parentArcs;
//	vector<uint> childArcs;
//	bool prune = false;
//	vector<uint> path;
//
//	DDNode(){} // ASAP complete this constructor.
//
//	DDNode(uint id): nodeId{id}{}
//	/*
//
//	DDNode(uint id, const forward_list<uint>& st, uint layerNumber): nodeId{id}, state{st},
//		layerNo{layerNumber}, childArcs{}, objectiveVal{0}, isExact{true}, parentId{0} {
//		if (std::find(st.begin(), st.end(), 0) == st.end()) state.push_front(0);
//	}
//
//	DDNode(uint32_t id, const forward_list<uint32_t>& st, uint16_t layerNumber, uint64_t objective, uint32_t parent, const vector<DDArc>& child)
//			: nodeId{id}, state{st}, layerNo{layerNumber}, objectiveVal{objective}, isExact{true}, parentId{parent}, childArcs{child} {
//		if (std::find(st.begin(), st.end(), 0) == st.end()) state.push_front(0);
//	}
//
//	DDNode(const DDNode& node) // copy constructor
//			: nodeId{node.nodeId}, state{node.state}, layerNo{node.layerNo},objectiveVal{node.objectiveVal},
//			isExact{node.isExact}, parentId{node.parentId}, childArcs{node.childArcs}, prune{node.prune} {}
//
//	 */
//	// operators
//
//	DDNode& operator=(const DDNode& node) {
//		this->layerNo = node.layerNo;
//		this->nodeId = node.nodeId;
//		this->objectiveVal = node.objectiveVal;
//		this->isExact = node.isExact;
//		this->state = node.state;
//		return *this;
//	}
//
//	DDNode& operator=(DDNode&& node) noexcept {
//		this->nodeId = node.nodeId;
//		this->layerNo = node.layerNo;
//		this->objectiveVal = node.objectiveVal;
//		this->isExact = node.isExact;
//		this->state = node.state;
//		return *this;
//	}
//
//	bool operator==(const DDNode& node) const{
//		return this->nodeId == node.nodeId; /* ASAP fix this*/
//	}
//
//	bool operator<(const DDNode& node){
//		return this->nodeId < node.nodeId; /* ASAP fix this */
//	}
//
//	bool operator<=(const DDNode& node){
//		return true; /* fix this ASAP */
//	}
//
//	bool operator>(const DDNode& node){
//		return true; /* fix this ASAP */
//	}
//
//	bool operator>=(const DDNode& node){
//		return true; /* fix this ASAP */
//	}
//
//	struct HashFunction {
//		size_t operator()(const DDNode& node) const {
//			return std::hash<int>()(static_cast<int>(node.nodeId));
//		}
//	};
//};
//
//
//
//typedef std::vector<uint> NodeLayer;
//typedef std::vector<NodeLayer> DDTree;
//
//// Class definition for Relaxed DD.
//class RelaxedDD{
//
//private:
//	uint32_t upperBound;
//	uint32_t maxWidth;
//	PruneStrategy strategy = TRAIL;
//	bool isExact;
//	NodeLayer cutset { };
//	DDTree tree;
//
//	/*
//	 * Merge nodes based on a selected criteria.
//	 */
//	void mergeNodes(Layer& layer) {
//		/* TODO; complete this function */
//		/* debug this function */
//	}
//
//	inline void pruneNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer);
//	//inline void buildNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer, uint32_t currentVariable);
//
//public:
//	inline void build(const Network& network, uint32_t startVariable = 0, uint32_t gamma = 9999); // todo: set gamma later.
//
//	NodeLayer getCutSet();
//};
//
//// class definition for Restricted DD
//class RestrictedDD{
//private:
//	uint64_t lowerBound;
//	PruneStrategy strategy = TRAIL;
//	DDTree tree = DDTree();
//	//void trimNodes(Layer& currentLayer);
//	inline void pruneNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer);
//	//void insertNode(Layer& currentLayer, DDNode&& node);
//	//static inline void buildNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer, uint32_t currentVariable);
//
//
//public:
//	bool isExact;
//	uint32_t maxWidth;
//	Layer cutset;
//
//	RestrictedDD() {
//	}
//
//	explicit RestrictedDD(uint32_t max_width): maxWidth{max_width}, isExact{true}, lowerBound{0} {
//
//	}
//
//	void build(const Network& network, uint32_t startVariable= 0, uint32_t gamma=99999);
//
//	uint32_t getLowerBound() {return this->lowerBound; }
//
//	Layer getCutSet() { return cutset; }
//
//	vector<uint32_t> getPath();
//
//	void refineDD(RestrictedDD& restrictedDD);
//};
//
//
//// Class definition for general DD. Intended as parent class for both
//// restricted DD class and relaxed DD class.
//class DD {
//
//private:
//    //DDNode root;
//	//vector<vector<DDNode>> tree;
//	DDTree tree;
//	std::map<uint, DDNode> nodes;
//	std::map<uint, DDArc> arcs;
//    uint64_t objective;
//	uint maxWidth = 0;
//
//public:
//
//	DD(){
//		//this->tree = vector<unordered_set<DDNode>>(coeff.size()+2);
//		// TODO: add child nodes to the root.
//	}
//
//    uint64_t getObjectiveVal() {
//        return this->objective;
//    }
//
//	void getPath() {
//		// this function computes the path from terminal node to the root.
//		// start from terminal node and pick a maximum value and select that node.
//		// vector of labels from terminal to root.
//		vector<uint32_t> path(this->tree.size());
//
//		for (auto it = this->tree.rbegin(); it != this->tree.rend(); ++it){
//
//		}
//
//		uint32_t maxZ = 0;
//		uint32_t  maxLength = 0;
//		uint32_t  parentId;
//		// find max length of path.
//		for (const auto& node: this->tree[this->tree.size()-2]){
//			if (node.childArcs[0].label > maxZ) {
//				maxZ = node.childArcs[0].label;
//				parentId = node.nodeId;
//			}
//		}
//
//		//  iterate through layers.
//		for (int layer = this->tree.size()-2; layer > 0; layer--){
//			const auto parentNode = this->tree[layer][parentId];
//
//			parentId = parentNode.parentId;
//		}
//	}
//
//	static void buildLayer(NodeLayer& currentLayer, NodeLayer& nextLayer, uint variable){
//		//
//
//	}
//	/**
//	 * Builds next layer in the tree from the current layer.
//	 */
//	static void buildNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer, uint32_t currentVariable){
//		// TODO: need to add currentVariable to every arc.
//
//		// ASAP definitely use actual Indexset instead of nodeIds, since '0' can also be a nodeId.
//		uint32_t nodeId = 0;
//		for (auto& node: currentLayer){
//			uint32_t arcId = 0;
//			for (auto i: node.state){
//				// create an arc
//				DDArc arc{arcId++, nodeId, i};
//				arc.weight = currentVariable;
//				auto newState = node.state;
//				if (i != 0) newState.remove(i); // remove only non-zero states.
//				if (std::find(newState.begin(), newState.end(), 0 ) == newState.end())
//					newState.push_front(0); // if '0' is not present in state, add it.
//				DDNode newNode {nodeId++,newState, node.layerNo+1};
//
//				node.childArcs.push_back(arc);
//				nextLayer.push_back(std::move(newNode));
//			}
//		}
//	}
//
//	/**
//	 * Builds the DD sub-tree from the current variable
//	 */
//	void buildSubTree(const Network& network, uint32_t startVariable = 0, uint32_t gamma = 999){
//		const auto Vbar = network.Vbar;
//		const auto networkNodes = network.networkNodes;
//		const auto networkArcs = network.networkArcs;
//
//		NodeLayer currentLayer;
//		bool first = true;
//
//		for (auto iterator = Vbar.begin()+startVariable; iterator != Vbar.end(); iterator++){
//			auto networkNode = networkNodes[*iterator];
//
//			forward_list<uint32_t> indexSet;
//			for (int i = networkNode.outNodeIds.size()-1; i >= 0; i--)
//				indexSet.push_front(networkNode.outNodeIds[i]+1);
//
//			indexSet.push_front(0);
//
//			if (currentLayer.empty() && first){ // build subtree for the starting node.
//				DDNode rootNode{0, indexSet, 0};
//				currentLayer.push_back(rootNode);
//				this->tree.push_back(currentLayer);
//				printLayer(currentLayer);
//				// iterate over incoming arcs and build layers.
//				for (auto incomingId: networkNode.inNodeIds){
//					NodeLayer nextLayer;
//					buildNextLayer(currentLayer, nextLayer, incomingId);
//					printLayer(nextLayer);
//					this->tree.push_back(nextLayer);
//					currentLayer = std::move(nextLayer);
//				}
//				first = false;
//			}
//			else { // second variable from start, update the status of nodes in current layer.
//				for(auto& node: currentLayer) node.state = indexSet;
//
//				// iterate over incoming arcs.
//				for (const auto incomingArc: networkNode.inNodeIds){
//					NodeLayer nextLayer;
//					buildNextLayer(currentLayer, nextLayer, incomingArc);
//					this->tree.push_back(nextLayer);
//					printLayer(nextLayer);
//					currentLayer = std::move(nextLayer);
//				}
//			}
//		}
//		// terminal layer.
//		NodeLayer terminalLayer;
//		DDNode terminalNode {0,{0},static_cast<uint32_t>(networkNodes.size()+1)};
//		terminalLayer.push_back(terminalNode);
//
//		for (auto& node: currentLayer){
//			DDArc arc1{0,0, gamma};
//			DDArc arc2{1,0,-gamma};
//			node.childArcs.push_back(arc1);
//			node.childArcs.push_back(arc2);
//		}
//		//printTree();
//	}
//	/**
//	 * Builds the DD tree given the network.
//	 */
//	void build(const Network& network, uint16_t startVariable=0, uint32_t gamma=999999){
//		// TODO set gamma.
//		// TODO: instead of starting from root, this function should also builds tree from the given layer.
//		const auto Vbar = network.Vbar;
//		const auto networkNodes = network.networkNodes;
//		const auto networkArcs = network.networkArcs;
//
//		NodeLayer currentLayer;
//
//		bool first = true;
//		for (const auto vBarId: Vbar){
//
//			auto networkNode = networkNodes[vBarId];
//			// index set of current node.
//			forward_list<uint32_t> indexSet;
//			for (int i = networkNode.outNodeIds.size() -1; i >= 0; i--)
//				indexSet.push_front(networkNode.outNodeIds[i]+1);
//			indexSet.push_front(0);
//
//			//for (const auto id: networkNode.outNodeIds)
//			//	indexSet.push_front(id+1); // indexset {nodeId_1+1, nodeId_2+1, ... }
//			//indexSet.push_front(0);
//			//forward_list<uint32_t> indexSet {networkNode.outNodeIds.begin(), networkNode.outNodeIds.end()};
//
//			if (currentLayer.empty() && first){
//				// build tree for first node in v_bar set the first to false.
//				// networkNode is root node. build a DDNode for root and push back to current layer.
//				DDNode rootNode{0, indexSet, 0};
//				currentLayer.push_back(rootNode);
//				this->tree.push_back(currentLayer); // first layer in tree.
//				// iterate over incoming arcs and build next layer.
//				for (auto incomingId: networkNode.inNodeIds){
//					NodeLayer nextLayer;
//					buildNextLayer(currentLayer, nextLayer, incomingId);
//					this->tree.push_back(nextLayer);
//					printLayer(nextLayer);
//					currentLayer = std::move(nextLayer);
//				}
//				first = false;
//			}
//			else { // update the state of nodes in current layer to state of new node in vbar.
//				for (auto& node: currentLayer){
//					node.state = indexSet;
//				}
//				// iterate over the incoming arcs of the current selected node in vbar.
//				for (const auto incomingArc: networkNode.inNodeIds){
//					NodeLayer nextLayer;
//					buildNextLayer(currentLayer, nextLayer, incomingArc);
//					this->tree.push_back(nextLayer);
//					printLayer(nextLayer);
//					currentLayer = std::move(nextLayer);
//				}
//			}
//		}
//		// terminal layer.
//		NodeLayer terminalLayer;
//		DDNode terminalNode{0,{0},100};
//		terminalLayer.push_back(terminalNode);
//
//		for (auto& node: currentLayer){
//			DDArc child1{0,0, gamma};
//			DDArc child2{1,0,-gamma};
//			node.childArcs.push_back(child1);
//			node.childArcs.push_back(child2);
//		}
//
//		printTree();
//	}
//
//	void buildOne(const Network& network) {
//
//		uint32_t label = 100;
//
//		const auto Vbar = network.Vbar;
//		const auto networkNodes = network.networkNodes;
//		const auto networkArcs = network.networkArcs;
//
//		// extract a node form vbar and set it as root.
//		auto id = Vbar[0];
//		auto networkNode = networkNodes[id];
//		forward_list<uint32_t> indexSet {networkNode.outNodeIds.begin(), networkNode.outNodeIds.end()};
//		// build a DDNode with it.
//		DDNode root{0,indexSet, 0};
//		// first layer
//		vector<DDNode> currentLayer;
//		currentLayer.push_back(root);
//
//		for (uint32_t i = 0; i < networkNode.inDegree; i++) {
//			vector<DDNode> nextLayer;
//			uint32_t idNewNode = 0;
//			// iterate over nodes in the current layer.
//			for (auto &node: currentLayer) {
//				// iterate through index set, create arcs and insert them to node.
//				uint32_t arcId = 0;
//				for (const auto &s: node.state) {
//					// create arc and a node
//					auto new_state = node.state;
//					new_state.remove(s);
//
//					DDArc arc{arcId++, idNewNode, s};
//					DDNode child{idNewNode++, new_state, node.layerNo+1};
//					node.childArcs.push_back(arc);
//					nextLayer.push_back(std::move(child));
//				}
//			}
//			tree.push_back(std::move(currentLayer));
//			currentLayer = nextLayer;
//		}
//		// for all the nodes in the last layer, add two arcs a_1, a_2 that connects terminal.
//		vector<DDNode> lastLayer;
//		DDNode terminal{0,{0},0};
//		lastLayer.push_back(terminal);
//		for (auto& node: currentLayer){
//			//
//			DDArc child1{0, 0, label};
//			DDArc child2{1, 0, -label};
//			node.childArcs.push_back(child1);
//			node.childArcs.push_back(child2);
//		}
//		printTree();
//	}
//
//	void printTree(){
//
//		for (const auto& layer: tree){
//			for (const auto& node: layer){
//				std::cout << node.nodeId << " ";
//			}
//			std::cout << std::endl;
//		}
//	}
//
//	void printLayer(const NodeLayer& layer){
//		for (const auto& node: layer){
//			std::cout << node.nodeId << " ";
//		}
//		std::cout << std::endl;
//	}
//};
//
//
//class DDClass {
//
//private:
//	DDType type;
//	DDTree tree;
//	std::unordered_map<uint, DDNode> nodes;
//	std::unordered_map<uint, DDArc> arcs;
//	uint maxWidth;
//	//DDNode root;
//	uint lastId = 0;
//	PruneStrategy strategy = PruneStrategy::TRAIL;
//
//	NodeLayer&& buildNextLayer(const NodeLayer& previousLayer, uint index);
//
//public:
//	DDClass(DDType type1): type{type1}, maxWidth{128} { // start from first variable.
//
//		// initialize all values
//
//	}
//
//	DDClass(DDType type1, const DDNode& rootNode): type{type1}{
//		maxWidth = 128;
//		lastId = rootNode.nodeId;
//		// update root node.
//
//	}
//
//	void build(const Network& network, DDNode& node /* starting node */, uint index);
//
//	void build(const Network& network); /* builds from starting node */
//
//
//};
//
////#endif //SGUFP_SOLVER_DD_H
//

