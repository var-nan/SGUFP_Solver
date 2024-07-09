//
// Created by nandgate on 6/3/24.
//

#include "DD.h"

// helper functions to display tree and layer.
static void printTree(const vector<vector<DDNode>>& tree){

	for (const auto& layer: tree){
		for (const auto& node: layer){
			std::cout << node.parentId << " ";
		}
		std::cout << std::endl;
	}
}

static void printLayer(const NodeLayer& layer){
	for (const auto& node: layer){
		std::cout << node.parentId << " ";
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
			newNode.parentId = node.nodeId;
			node.childArcs.push_back(arc);
			nextLayer.push_back(std::move(newNode));
		}
	}
}


void buildTree(const Network& network, uint32_t startVariable, uint32_t gamma,
                      const std::function<void(NodeLayer&, NodeLayer&)>& pruneFunction, vector<vector<DDNode>> tree) {

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
			tree.push_back(currentLayer);
			//printLayer(currentLayer);
			// iterate over incoming arcs and build layers.
			for (auto incomingId: networkNode.inNodeIds) {
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingId);
				pruneFunction(currentLayer, nextLayer); // prune next layer.
				printLayer(nextLayer);
				tree.push_back(nextLayer);
				currentLayer = std::move(nextLayer);
			}
			first = false;
		} else { // second variable from start, update the status of nodes in current layer.
			for (auto &node: currentLayer) node.state = indexSet;

			// iterate over incoming arcs.
			for (const auto incomingArc: networkNode.inNodeIds) {
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingArc);
				pruneFunction(currentLayer, nextLayer);
				tree.push_back(nextLayer);
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
	tree.push_back(terminalLayer);
}



/**
 * Builds Restricted subtree starting the given variable as root.
 * Can also be called to build subtree from any variable in the defined order.
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
		for (int i = networkNode.outDegree - 1; i >= 0; i--) // insert elements from last.
			indexSet.push_front(networkNode.outNodeIds[i] + 1);
		indexSet.push_front(0); // add zero to index set.

		if (currentLayer.empty() && first) { // build subtree for the starting node.
			DDNode rootNode{0, indexSet, 0};
			currentLayer.push_back(rootNode);
			this->tree.push_back(currentLayer);
			//printLayer(currentLayer);
			// iterate over incoming arcs and build layers.
			for (auto incomingId: networkNode.inNodeIds) { // actual Id of the node.
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingId);
				if (nextLayer.size() > maxWidth)
					pruneNextLayer(currentLayer, nextLayer); // prune next layer.
				//printLayer(nextLayer);
				this->tree.push_back(nextLayer);
				currentLayer = std::move(nextLayer);
			}
			first = false;
		} else { // second variable from start, update the states of nodes in current layer to index set.
			for (auto &node: currentLayer) node.state = indexSet;

			// iterate over incoming arcs.
			for (const auto incomingArc: networkNode.inNodeIds) {
				NodeLayer nextLayer;
				buildNextLayer(currentLayer, nextLayer, incomingArc);
				if (nextLayer.size() > maxWidth) // prune next layer if maxWidth is reached.
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
		//DDArc arc2{1, 0, -gamma}; // I don't think this arc is necessary here.
		node.childArcs.push_back(arc1);
		//node.childArcs.push_back(arc2);
	}
	this->tree.push_back(terminalLayer);

	std::cout << this->tree.size() << std::endl;
	//RestrictedDD resDd;
	//auto fp = std::bind(&RestrictedDD::pruneNextLayer, *this, placeholders::_1, placeholders::_2);
	//buildTree(network, startVariable, gamma, fp, tree);

	//printTree(this->tree);
}

vector<uint32_t> RestrictedDD::getPath() {
	// return a path from root to terminal.
	// start from terminal node and pick a maximum value in terminal layer.
	vector<uint32_t> path(this->tree.size() -1);

	uint32_t maxZ = 0;
	uint32_t maxLength = 0;
	uint32_t parentId;

	auto terminalIndex = this->tree.size()-2;

	for (const auto& node: this->tree[terminalIndex]){
		if (node.childArcs[0].label > maxZ){
			maxZ = node.childArcs[0].label;
			parentId = node.nodeId;
		}
	}

	path[terminalIndex] = parentId;

	// iterate through other layers.
	for (auto layer = this->tree.size()-2; layer > 0; layer--){
		const auto parentNode = tree[layer][parentId];
		parentId = parentNode.parentId;
		path[layer] = parentNode.parentId;
	}

	return path;
}

DDNode& getNodeWithIndex(vector<DDNode>& layer, uint32_t index, bool& valid){
	// returns the node that has the given index as nodeId.
	for (auto& node: layer){
		if (node.nodeId == index) {
			valid = true;
			return node;
		}
	}
	valid = false;
}

bool isArcPresent(const vector<DDArc>& childArcs, const DDArc& arc){
	// return true if arc is present in the childArcs.
	for (auto& a: childArcs){
		if (a.label == arc.label)
			return true;
	}
	return false;
}

void RestrictedDD::refineDD(RestrictedDD &restrictedDD) {
	// start with the root node and node the process each layer.

	for (unsigned int i = 0; i < this->tree.size(); i++) {
		// start from the root thisTreeLayer, and go through all the
		auto& thisTreeLayer = this->tree[i];
		auto& otherTreeLayer = restrictedDD.tree[i];
		vector<bool> filter(thisTreeLayer.size(), false);
		unsigned int index = 0;

		for (auto& node: thisTreeLayer){
			bool valid = false;
			auto& otherNode = getNodeWithIndex(otherTreeLayer, node.nodeId, valid);
			if (valid){ // node is present in both layers.
				/*
				for (auto& arc: otherNode.childArcs){ // TODO: delete arcs on the fly?

					// i love you, its ruining my life, i touch you for only a fortnight.
					// if arc is not present in the other node, then remove it.
					// iterate through childArcs and remove
					auto it = find_if(node.childArcs.begin(), node.childArcs.end(), [&arc](const DDArc& a){return a.label == arc.label;});
					if (it == node.childArcs.end()){
						// remove that arc.
					}
				}
				*/
				// intersection of arcs of both nodes.
				vector<DDArc> tempArcs;
				for (const auto& arc: node.childArcs){
					if (isArcPresent(otherNode.childArcs, arc))
						tempArcs.push_back(arc);
				}
				node.childArcs = tempArcs;
				filter[index] = true;
			}
			index++;
		}

		// remove nodes in the layer.
		// ASAP remove unfiltered nodes in the current layer.
		// update the child arcs of parent nodes.
		// update the index of remaining child nodes for its parents.

	}
}

// INFO: RELAXED DECISION DIAGRAM.

void RelaxedDD::pruneNextLayer(NodeLayer &currentLayer, NodeLayer &nextLayer) {

	if (strategy == TRAIL){
		// ASAP check maxWidth index.
		// merge all trailing nodes to single node (last node in nextLayer).
		unordered_set<uint32_t> states;

		for (int i = maxWidth-1; i < nextLayer.size(); i++){
			auto state = nextLayer[i].state;
			states.insert(state.begin(), state.end());
		}

		forward_list<uint32_t> newState {states.begin(), states.end()};

		DDNode newNode{maxWidth-1, newState, nextLayer[0].layerNo}; //

		NodeLayer newNextLayer(nextLayer.begin(), nextLayer.begin()+(maxWidth-2));
		newNextLayer.push_back(newNode);
		nextLayer = std::move(newNextLayer);

		// update arcs in current layer.
		uint32_t i = 0;
		for (auto& node: currentLayer){
			for (auto& arc: node.childArcs){
				if (i >= maxWidth) // all arcs point to last node in nextLayer.
					arc.headId = maxWidth-1; // TODO: check if is maxWidth-1 or maxWidth-2;
				i++;
			}
		}
	}
}

void RelaxedDD::build(const Network &network, uint32_t startVariable, uint32_t gamma) {

}