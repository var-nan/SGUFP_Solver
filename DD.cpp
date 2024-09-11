//
// Created by nandgate on 6/3/24.
//


#include "DD.h"

void DD::build(const vector<pair<int,int>> processingOrder, DDNode &node, int index) {

	// set the node to the root and start from there.

	vector<int> currentLayer; // should be base layer.
	currentLayer.push_back(node.id);

	//int lastInserted = 0;
	//tree.push_back({0});

	DDNode root(node); // need to reset the node Id to zero or 1.
	root.outgoingArcs.clear();
	root.incomingArcs.clear();
	// insert root to nodes.
	tree.push_back(currentLayer);
	auto start = processingOrder.begin()+index;
	auto end = processingOrder.end();

	for (; start < end; start++){

		auto[a,b] = *start;
		vector<int>&& nextLayer = buildNextLayer(currentLayer, index);
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}
}

vector<int>&& DD::buildNextLayer(vector<int>& currentLayer, int index) {
	/*
	 * build next layer from the given current layer.
	 * should also create arcs and store them in the map.
	 * should also create nodes and store them in the map.
	 */

	vector<int> nextLayer;
	//set<pair<set<int>,int>> allStates;

	for (auto id: currentLayer){ // iterate through current layer.

		auto& node = nodes[id];

		// add zero state
		DDNode zeroNode;
		zeroNode.id = lastNode++; // TODO: need to add some more fields.
		zeroNode.states = node.states;
		DDArc zeroArc(lastArc++, zeroNode.id, id, 0);
		node.outgoingArcs.push_back(zeroArc.id);
		zeroNode.incomingArcs.push_back(zeroArc.id);
		// insert to map.
		nodes.insert({zeroNode.id, zeroNode});
		arcs.insert({zeroArc.id, zeroArc});

		nextLayer.push_back(zeroNode.id);

		for (auto state: node.states){
			// create arc and node
			DDNode node1;// create a new node with some Id
			node1.id = lastNode++;
			auto states = node.states; // remove the state from states.
			states.erase(std::remove(states.begin(), states.end(), state), states.end());
			node1.states = std::move(states);
			DDArc arc(lastArc++, node1.id, id, state);
			node1.incomingArcs.push_back(arc.id);
			node.outgoingArcs.push_back(arc.id); // NOTE: SHOULD BE ID OF THE CHILDARC, NOT CHILD NODE.
			//DDArc arc(lastArc++, 10, 10, 10);
			nodes.insert({node1.id, node1});
			arcs.insert({arc.id, arc});
			nextLayer.push_back(node1.id); // ASAP: remove the state from the list of states and
		}
	}

	return std::move(nextLayer);
}

void DD::mergeNodes(DDNode& node1, DDNode& node2) {
	/*
	 * Merge node2 with node1. Updates node1 attributes and removes node2 from dictionary.
	 * Used in reduceNodes().
	 */
	for (auto parent: node2.incomingArcs){
		auto& arc = arcs[parent];
		arc.head = node1.id; // TODO: change other attributes if needed
	}

	for (auto child: node2.outgoingArcs){
		auto& arc = arcs[child];
		arc.tail = node1.id; // TODO: change other attributes if needed.
	}
	// update node1 incomingArcs to add node2's incomingArcs.
	node1.incomingArcs.insert(node1.incomingArcs.end(), node2.incomingArcs.begin(), node2.incomingArcs.end());
	node1.outgoingArcs.insert(node1.outgoingArcs.end(), node2.outgoingArcs.begin(), node2.outgoingArcs.end());

	node2.incomingArcs.clear();
	node2.outgoingArcs.clear();

}

void DD::reduceLayer(vector<uint>& currentLayer) {
	/*
	 * Merge two nodes that containing same state. Update incoming and outgoing arcs to point to merged node.
	 */

	vector<bool> filter(currentLayer.size()); // to keep track of current layer nodes.
	unordered_set<tuple<unordered_set<int>,int,int>,tuple_hash,tuple_equal> allStates;
	uint j = 0;

	for (auto i: currentLayer){
		auto& node = nodes[i];
		auto [it,isInserted] = allStates.insert({node.states, node.state2, i});

		if (!isInserted) { // already exists, merge both nodes.
			auto& [temp_state, temp_stat2, pos ] = *(it);
			// get node from dict with pos as key and update its attributes.
			auto& node2 = nodes[pos];
			mergeNodes(node2, node);
			// delete node from nodes.
			nodes.erase(node.id);
			filter[j] = true;
		}
		j++;
	}
	// remove filtered nodeIds from current layer.
	currentLayer.erase(std::remove_if(currentLayer.begin(), currentLayer.end(),
	                                  [&filter](int x){return filter[x];}),
	                   currentLayer.end());

}

void DD::deleteArcById(uint id){
	/*
	 * deletes arc by id.`
	 * removes the arc from tail and head nodes.
	 * deletes arc from arcs dictionary.
	 *
	 * INFO function assumes both the head and tail nodes are in the nodes dictionary.
	 */
	auto& arc = arcs[id];
	auto& tailNode = nodes[arc.tail];
	auto& headNode = nodes[arc.head];
	headNode.incomingArcs.erase(std::find(headNode.incomingArcs.begin(), headNode.incomingArcs.end(), id));
	tailNode.outgoingArcs.erase(std::find(tailNode.outgoingArcs.begin(), tailNode.outgoingArcs.end(), id));
	arcs.erase(id);
}

void DD::deleteNodeById(uint id) {
	/*
	 * deletes node by its id.
	 * delete outgoing arcs vector. destructor deletes remaining attributes.
	 */
	auto& node = nodes[id];
	auto outgoingArcs = node.outgoingArcs;
	for (const auto arcId : outgoingArcs){
		deleteArcById(arcId);
	}
}

inline DDNode duplicate(const DDNode& node){
	// clone the node. copy outgoing arcs. do not copy incoming arcs.

}

void DD::duplicateNode(uint id){
	/*
	 * for every incoming arc, create new node (copy state parameters)
	 *  TODO: complete it today.
	 */
	auto& node = nodes[id];

	auto incomingArcs = node.incomingArcs;

	for (uint i = 1; i < incomingArcs.size(); i++) {
		//
		int arcId = incomingArcs[i];
		DDArc& arc = arcs[arcId];

		// TODO compute state.
		bool feasible = true; // TODO: update it.
		if (feasible){
			DDNode newNode = duplicate(node); // set incoming arc.
			newNode.incomingArcs.push_back(arcId);
			nodes.insert({newNode.id, newNode});
			// update head of arc to point to new node.
			arc.head = newNode.id;
		}
		else {
			// delete arc.
			deleteArcById(arcId);
			auto& tailOutArcs = nodes[arc.tail].outgoingArcs;
			tailOutArcs.erase(std::find(tailOutArcs.begin(), tailOutArcs.end(), arcId));
			// remvoe arc
			arcs.erase(arcId);
		}
	}
	// remove the incoming arcs of original node.
	// compute state
	if (true){ // feasible, delete all incoming arcs except the first.
		node.incomingArcs.erase(node.incomingArcs.begin()+1, node.incomingArcs.end());
	}
	else {
		// delete node.
		deleteArcById(node.incomingArcs[0]);
		deleteNodeById(node.id);
	}
}

static void printTree(const vector<vector<int>>& tree){

	for (const auto& layer: tree){
		for (const int node: layer){
			cout << node << " ";
		}
		cout << endl;
	}
}

int main(){
	// some processing order
	vector<pair<int,int>> processingOrder = {{1,2},{2,3}};

	DD diagram;

	DDNode rootNode;
	rootNode.states = unordered_set<int>({2,3,4});

	diagram.build(processingOrder, rootNode, 0);

	printTree(diagram.tree);
}
