//
// Created by nandgate on 6/3/24.
//

#include "DD.h"

// helper functions to display tree and layer.
static void printTree(const vector<vector<DDNode>>& tree){

	for (const auto& layer: tree){
		for (const auto& node: layer){
			std::cout << node.nodeId << " ";
		}
		std::cout << std::endl;
	}
}

static void printLayer(const NodeLayer& layer){
	for (const auto& node: layer){
		std::cout << node.nodeId << " ";
	}
	std::cout << std::endl;
}


/**
 * Prunes the given layers by maxWidth parameter.
 * Nodes are removed in the 'nextLayer' object and their corresponding Arcs are
 * updated in the nodes of "currentLayer" object.
 */
inline void RestrictedDD::pruneNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer){

	if (strategy == TRAIL){
		// iterate through all arcs in the current layer with a counter and remove all the arcs past the maxWidth.
		uint32_t i = 0;
		for (auto& node: currentLayer) {
			if (i >= maxWidth) node.childArcs.clear(); // remove all the child arcs of node.
			else {
				for (uint32_t arcIndex = 0; arcIndex < node.childArcs.size(); arcIndex++) {
					if (i >= maxWidth) {
						node.childArcs.erase(node.childArcs.begin() + arcIndex, node.childArcs.end());
						break;
					}
					i++;
				}
			}
			/*
			for (auto& arc: node.childArcs){
				if (i > maxWidth){
					//node.childArcs.erase(std::remove(std::execution::seq, std::begin(node.childArcs), std::end(node.childArcs), arc) , node.childArcs.end());
					//(std::find(node.childArcs.begin(), node.childArcs.end(), arc) != std::end(node.childArcs));
				}
			}
			 */
		}
		// remove all the trailing nodes in the next layer.
		nextLayer.erase(nextLayer.begin()+maxWidth, nextLayer.end());
	}
	else if (strategy == RANDOM) {
		// LATER:
	}
}

void NextLayer(NodeLayer &currentLayer, NodeLayer& nextLayer) {
	// iterate through all nodes and build next layer
	//NodeLayer nextLayer;
	uint32_t nodeId = 0;

	for (auto& node: currentLayer){
		uint32_t arcId = 0;
		for (const auto s: node.state) {
			DDArc arc{arcId++, nodeId, s};
			auto newState = node.state; // TODO: remove the current selected attribute and add '0'.
			newState.remove(s);
			DDNode newNode{nodeId++, newState, node.layerNo+1};
			node.childArcs.push_back(arc);
			nextLayer.push_back(std::move(newNode));
		}
	}

	//pruneNextLayer(currentLayer, nextLayer);

	//return std::move(nextLayer);
}

/**
 * Builds next layer with the current layer.
 */
static void buildNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer, uint32_t currentVariable){

	uint32_t nodeId = 0;
	for (auto& node: currentLayer){
		uint32_t arcId = 0;
		for (auto i: node.state){
			// create an arc and add it to child nodes.
			DDArc arc{arcId++, nodeId, i};
			arc.weight = currentVariable;
			auto newState = node.state;
			if (i != 0) newState.remove(i); // remove only non-zero states.
			//if (std::find(newState.begin(), newState.end(), 0 ) == newState.end())
			//	newState.push_front(0); // if '0' is not present in state, add it.
			DDNode newNode {nodeId++,newState, node.layerNo+1};
			node.childArcs.push_back(arc);
			nextLayer.push_back(std::move(newNode));
		}
	}
}


/**
 * Builds Restricted subtree starting the given variable as root.
 */
void RestrictedDD::build(const Network& network, uint32_t startVariable, uint32_t gamma) {
	const auto Vbar = network.Vbar;
	const auto networkNodes = network.networkNodes;
	const auto networkArcs = network.networkArcs;

	NodeLayer currentLayer;

	bool first = true;

	for (auto iterator = Vbar.begin() + startVariable; iterator != Vbar.end(); iterator++) {
		auto networkNode = networkNodes[*iterator];

		forward_list<uint32_t> indexSet;
		for (int i = networkNode.outDegree - 1; i >= 0; i--)
			indexSet.push_front(networkNode.outNodeIds[i] + 1);
		indexSet.push_front(0);

		if (currentLayer.empty() && first) { // build subtree for the starting node.
			DDNode rootNode{0, indexSet, 0};
			currentLayer.push_back(rootNode);
			this->tree.push_back(currentLayer);
			//printLayer(currentLayer);
			// iterate over incoming arcs and build layers.
			for (auto incomingId: networkNode.inNodeIds) {
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingId);
				pruneNextLayer(currentLayer, nextLayer); // prune next layer.
				//printLayer(nextLayer);
				this->tree.push_back(nextLayer);
				currentLayer = std::move(nextLayer);
			}
			first = false;
		} else { // second variable from start, update the status of nodes in current layer.
			for (auto &node: currentLayer) node.state = indexSet;

			// iterate over incoming arcs.
			for (const auto incomingArc: networkNode.inNodeIds) {
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingArc);
				pruneNextLayer(currentLayer, nextLayer);
				this->tree.push_back(nextLayer);
				printLayer(nextLayer);
				currentLayer = std::move(nextLayer);
			}
		}
	}
	// terminal layer.
	NodeLayer terminalLayer;
	DDNode terminalNode{0, {0}, static_cast<uint32_t>(networkNodes.size() + 1)};
	terminalLayer.push_back(terminalNode);
	// add arcs to terminal node.
	for (auto &node: currentLayer) {
		DDArc arc1{0, 0, gamma};
		DDArc arc2{1, 0, -gamma};
		node.childArcs.push_back(arc1);
		node.childArcs.push_back(arc2);
	}
}


void RelaxedDD::pruneNextLayer(NodeLayer &currentLayer, NodeLayer &nextLayer) {

	if (strategy == TRAIL){
		// merge all trailing nodes to single node.
		uint32_t i = 0;

		for (auto& node: currentLayer){
			// if iteration reached maxWidth, remove all the trailing arcs.
			if (i >= maxWidth) node.childArcs.clear();
			else {
				for (uint32_t arcIndex = 0; arcIndex < node.childArcs.size(); arcIndex++){
					if (i >= maxWidth){
						node.childArcs.erase(node.childArcs.begin()+arcIndex, node.childArcs.end());
						break;
					}
					i++;
				}
			}
		}

		// update next layer nodes.
		unordered_set<uint32_t> states;

		for (i = maxWidth; i < nextLayer.size(); i++){
			auto state = nextLayer[i].state;
			states.insert(state.begin(), state.end());
		}

		forward_list<uint32_t> newState {states.begin(), states.end()};
		DDNode newNode{maxWidth, newState, nextLayer[0].layerNo};

		NodeLayer newNextLayer(nextLayer.begin(), nextLayer.begin()+maxWidth);
		newNextLayer.push_back(newNode);
	}
}
