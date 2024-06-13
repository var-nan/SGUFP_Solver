//
// Created by nandgate on 6/1/24.
//

#ifndef SGUFP_SOLVER_DD_H
#define SGUFP_SOLVER_DD_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <unordered_set>
#include "Network.h"

using namespace std;

typedef vector<uint32_t> vui;
typedef vector<int> vi;

class DDArc{
public:
	uint32_t id;
	uint32_t headId; // id of child node in the next layer of tree.
	uint32_t label; // decision.
	int weight;

	DDArc(uint32_t arc_id, uint32_t tail_id, uint32_t arc_label)
		: id{arc_id}, headId{tail_id}, label{arc_label} {}
};

class DDNode {
public:
	uint32_t nodeId; /* id of the node */
	uint64_t objectiveVal; /* objective value is an unsigned long number */
	forward_list<uint32_t> state; /* node might have multiple state values when refining */
	uint32_t layerNo; /* with respect to global ordering of the variables */
	bool isExact;
	uint32_t parentId;
	vector<DDArc> childArcs;

	DDNode(){}; // ASAP complete this constructor.

	DDNode(uint32_t id, forward_list<uint32_t> st, uint32_t layerNumber): nodeId{id}, state{st}, layerNo{layerNumber}, childArcs{} {}

	DDNode(forward_list<uint32_t> st, uint16_t layerNumber, uint64_t objective, uint32_t parent, vector<DDArc> child)
			: state{st}, layerNo{layerNumber}, objectiveVal{objective}, isExact{true}, parentId{parent}, childArcs{child} {}

	DDNode(const DDNode& node) // TODO: copy parent and children to the new node.
			: state{node.state}, layerNo{node.layerNo}, objectiveVal{node.objectiveVal}, isExact{node.isExact} {}

	/*
	 * operators
	 */

	DDNode& operator=(const DDNode& node) {

		return *this;
	}

	DDNode& operator=(DDNode&& node) noexcept {
		this->layerNo = node.layerNo;
		this->objectiveVal = node.objectiveVal;
		this->isExact = node.isExact;
		this->state = node.state;
		return *this;
	}

	bool operator==(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator<(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator<=(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator>(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator>=(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	struct HashFunction {
		size_t operator()(const DDNode& node) const {
			return std::hash<int>()(static_cast<int>(node.nodeId));
		}
	};
};


typedef std::unordered_set<DDNode, DDNode::HashFunction> Layer;

class RelaxedDD{

private:
	uint32_t upperBound;
	uint32_t maxWidth;
	bool isExact;
	Layer cutset { };
	
	/*
	 * Merge nodes based on a selected criteria.
	 */
	void mergeNodes(Layer& layer) {
		/* TODO; complete this function */
		/* debug this function */
	}
	
	void insertNode(Layer& currentLayer, DDNode&& node){
	
	}
	
public:


};

class RestrictedDD{
private:
	uint64_t lowerBound;
	void trimNodes(Layer& currentLayer);

	void insertNode(Layer& currentLayer, DDNode&& node);

public:
	bool isExact;
	uint32_t maxWidth;
	Layer cutset;

	RestrictedDD(uint32_t max_width): maxWidth{max_width}, isExact{true}, lowerBound{0} {}

	void build(const DDNode& root, const vui& coefficients, const vui& values);

	uint32_t getLowerBound() {return this->lowerBound; }

	Layer getCutSet() { return cutset; }
};

typedef std::vector<DDNode> NodeLayer;

class DD {

private:
    DDNode root;
	vector<vector<DDNode>> tree;
    uint64_t objective;

public:

	DD(): tree{10}{
		//this->tree = vector<unordered_set<DDNode>>(coeff.size()+2);
		// TODO: add child nodes to the root.
	}

    uint64_t getObjectiveVal() {
        return this->objective;
    }

	void getPath() {
		// this function computes the path from terminal node to the root.
	}

	static NodeLayer&& buildNextLayer(NodeLayer& currentLayer, uint32_t currentVariable){
		// TODO: need to add currentVariable to every arc.

		// takes current layer and builds next layer based on the indexset.
		NodeLayer nextLayer;
		uint32_t nodeId = 0;

		for (auto& node: currentLayer){
			uint32_t arcId = 0;
			for (auto i: node.state){
				// create an arc
				DDArc arc{arcId++, nodeId, i};
				auto newState = node.state;
				newState.remove(i); // TODO: handle '0' label.
				DDNode newNode {nodeId++,newState, node.layerNo+1};
				node.childArcs.push_back(arc);
				nextLayer.push_back(std::move(newNode));
			}
		}
		return std::move(nextLayer);
	}

	void build(const Network& network, uint16_t startLayer=0, double gamma=-999999){
		// TODO: set gamma later.
		const auto Vbar = network.Vbar;
		const auto networkNodes = network.networkNodes;
		const auto networkArcs = network.networkArcs;

		NodeLayer currentLayer;

		bool first = true;

		for (const auto vBarId: Vbar){
			auto networkNode = networkNodes[vBarId];
			// index set of current node.
			forward_list<uint32_t> indexSet {networkNode.outArcIds.begin(), networkNode.outArcIds.end()};

			if (currentLayer.empty() && first){
				// build tree for first node in v_bar set the first to false.
				// networkNode is root node. build a DDNode for root and push back to current layer.
				DDNode rootNode{0, indexSet, 0};
				currentLayer.push_back(rootNode);
				// iterate over incoming arcs and build next layer
				for (auto incomingId: networkNode.inArcIds){
					auto nextLayer = buildNextLayer(currentLayer, incomingId);
					this->tree.push_back(nextLayer);
					currentLayer = std::move(nextLayer);
				}
				first = false;
			}
			else { // update the state of nodes in current layer to state of new node in vbar.
				for (auto& node: currentLayer){
					node.state = indexSet;
				}
				// iterate over the incoming arcs of the current selected node in vbar.
				for (const auto incomingArc: networkNode.inArcIds){
					auto nextLayer = buildNextLayer(currentLayer, incomingArc); // TODO: don't need to pass index set here, since nodes in current layer already have new state.
					this->tree.push_back(nextLayer);
					currentLayer = std::move(nextLayer);
				}
			}
		}
		// TODO: add terminal node.
		NodeLayer terminalLayer;
		DDNode terminalNode;
		terminalLayer.push_back(terminalNode);

		for (auto& node: currentLayer){
			DDArc child1{0,0, gamma};
			DDArc child2{1,0,-gamma};
			node.childArcs.push_back(child1);
			node.childArcs.push_back(child2);
		}
	}

	void buildOne(const Network& network) {

		uint32_t label = 100;

		const auto Vbar = network.Vbar;
		const auto networkNodes = network.networkNodes;
		const auto networkArcs = network.networkArcs;

		// extract a node form vbar and set it as root.
		auto id = Vbar[0];
		auto networkNode = networkNodes[id];
		forward_list<uint32_t> indexSet {networkNode.outArcIds.begin(), networkNode.outArcIds.end()};
		// build a DDNode with it.
		DDNode root{0,indexSet, 0};
		// first layer
		vector<DDNode> currentLayer;
		currentLayer.push_back(root);

		for (uint32_t i = 0; i < networkNode.inDegree; i++) {
			vector<DDNode> nextLayer;
			uint32_t idNewNode = 0;
			// iterate over nodes in the current layer.
			for (auto &node: currentLayer) {
				// iterate through index set, create arcs and insert them to node.
				uint32_t arcId = 0;
				for (const auto &s: node.state) {
					// create arc and a node
					auto new_state = node.state;
					new_state.remove(s);

					DDArc arc{arcId++, idNewNode, s};
					DDNode child{idNewNode++, new_state, node.layerNo+1};
					node.childArcs.push_back(arc);
					nextLayer.push_back(std::move(child));
				}
			}
			tree.push_back(std::move(currentLayer));
			currentLayer = nextLayer;
		}
		// for all the nodes in the last layer, add two arcs a_1, a_2 that connects terminal.
		vector<DDNode> lastLayer;
		DDNode terminal{0,0,0};
		lastLayer.push_back(terminal);
		for (auto& node: currentLayer){
			//
			DDArc child1{0, 0, label};
			DDArc child2{1, 0, -label};
			node.childArcs.push_back(child1);
			node.childArcs.push_back(child2);
		}
	}
};



#endif //SGUFP_SOLVER_DD_H
