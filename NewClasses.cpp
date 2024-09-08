//
// Created by nandgate on 9/5/24.
//

#include "NewClasses.h"

void Diagram::build(const vector<pair<int,int>> processingOrder, Node &node, int index) {

	// set the node to the root and start from there.

	vector<int> currentLayer; // should be base layer.
	currentLayer.push_back(node.id);

	//int lastInserted = 0;
	//tree.push_back({0});

	Node root(node);
	root.children.clear();
	root.parents.clear();
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

vector<int>&& Diagram::buildNextLayer(vector<int>& currentLayer, int index) {
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
		Node zeroNode;
		zeroNode.id = lastNode++; // TODO: need to add some more fields.
		zeroNode.states = node.states;
		Arc zeroArc(lastArc++, zeroNode.id, id, 0);
		node.children.push_back(zeroArc.id);
		zeroNode.parents.push_back(zeroArc.id);
		// insert to map.
		nodes.insert({zeroNode.id, zeroNode});
		arcs.insert({zeroArc.id, zeroArc});

		nextLayer.push_back(zeroNode.id);

		for (auto state: node.states){
			// create arc and node
			Node node1;// create a new node with some Id
			node1.id = lastNode++;
			auto states = node.states; // remove the state from states.
			states.erase(std::remove(states.begin(), states.end(), state), states.end());
			node1.states = std::move(states);
			Arc arc( lastArc++, node1.id, id, state);
			node1.parents.push_back(arc.id);
			node.children.push_back(arc.id); // NOTE: SHOULD BE ID OF THE CHILDARC, NOT CHILD NODE.
			//Arc arc(lastArc++, 10, 10, 10);
			nodes.insert({node1.id, node1});
			arcs.insert({arc.id, arc});
			nextLayer.push_back(node1.id); // ASAP: remove the state from the list of states and
		}
	}

	return std::move(nextLayer);
}

void Diagram::mergeNodes(Node& node1, Node& node2) {
	/*
	 * Merge node2 with node1. Updates node1 attributes and removes node2 from dictionary.
	 * Used in reduceNodes().
	 */
	for (auto parent: node2.parents){
		auto& arc = arcs[parent];
		arc.head = node1.id; // TODO: change other attributes if needed
	}

	for (auto child: node2.children){
		auto& arc = arcs[child];
		arc.tail = node1.id; // TODO: change other attributes if needed.
	}
	// update node1 parents to add node2's parents.
	node1.parents.insert(node1.parents.end(), node2.parents.begin(), node2.parents.end());
	node1.children.insert(node1.children.end(), node2.children.begin(), node2.children.end());

	node2.parents.clear();
	node2.children.clear();

}

void Diagram::reduceLayer(vector<uint>& currentLayer) {
	/*
	 * Merge two nodes that containing same state. Update incoming and outgoing arcs to point to merged node.
	 */

	vector<bool> filter(currentLayer.size());
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

	Diagram diagram;

	Node rootNode;
	rootNode.states = unordered_set<int>({2,3,4});

	diagram.build(processingOrder, rootNode, 0);

	printTree(diagram.tree);
}