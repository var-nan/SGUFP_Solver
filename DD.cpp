//
// Created by nandgate on 6/3/24.
//

#include "DD.h"

void RestrictedDD::pruneNextLayer(NodeLayer& currentLayer, NodeLayer& nextLayer){

	if (strategy == TRAIL){
		// remove all the trailing nodes in the next layer.
		nextLayer.erase(nextLayer.begin()+maxWidth, nextLayer.end());
		// update arcs in the current layer.
		uint32_t i = 0;
		for (auto& node: currentLayer) {
			for (auto& arc: node.childArcs){
				if (i > 32){
					//node.childArcs.erase(std::remove(std::execution::seq, std::begin(node.childArcs), std::end(node.childArcs), arc) , node.childArcs.end());
					//(std::find(node.childArcs.begin(), node.childArcs.end(), arc) != std::end(node.childArcs));
				}
			}
		}
	}
	else if (strategy == RANDOM) {
		// if random then don't know
	}
}

NodeLayer&& RestrictedDD::buildNextLayer(NodeLayer &currentLayer) {
	// iterate through all nodes and build next layer
	NodeLayer nextLayer;
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

	pruneNextLayer(currentLayer, nextLayer);

	return std::move(nextLayer);
}

void RestrictedDD::build(const DDNode& root, const vui& coefficients, const vui& values){

	//Layer currentLayer{root};

	int N = coefficients.size();
	
}