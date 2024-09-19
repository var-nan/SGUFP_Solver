//
// Created by nandgate on 6/3/24.
//


#include "DD.h"

/*
 * 1. change build next layer function to make the child not to store the solution vector.
 * 2. change build function to account for cutset generation.
 */


void DD::build(const Network& network, DDNode& node, int index) {
	// node parameter should be initialized before calling this function. node should contain states.
	// this node will be inserted as the root of the tree.

	// should reset nodes and arcs maps.
	const auto& arcOrder = network.processingOrder;
	const auto& stateUpdateMap = network.stateUpdateMap;

	// set the node to the root and start from there.
	vector<int> currentLayer; // should be base layer.

	startTree = index; // TODO; verify next.

	node.incomingArcs.clear();
	node.outgoingArcs.clear();
	node.id = 0; // make new node if necessary.
	nodes.insert(std::make_pair(node.id, node));
	currentLayer.push_back(node.id);
	// insert root layer to tree.
	tree.push_back(currentLayer);

	int i = 1;

	auto start = arcOrder.begin() + index;
	auto end = arcOrder.end();

	for (; start < end; start++){

		auto[a,b] = *start;
		if (stateUpdateMap.count(a)) { // update state of each node in the layer.
			const auto& newStates = stateUpdateMap.at(a);
			updateState(currentLayer, newStates);
		}
		//const unordered_set<int> temp = stateUpdateMap.at(a);
		i++;
		vector<int> nextLayer;
		bool isExact = buildNextLayer(currentLayer, nextLayer, index);
		if (isExact) exactLayer = i; // at last, this number should be exact layer number.
		//reduceLayer(nextLayer); // INFO not doing reduction.
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}
	// TODO: add terminal node layer.
}

inline void DD::updateState(const vector<int> &currentLayer, const unordered_set<int> &states){

	for (auto id: currentLayer){
		auto& node = nodes[id];
		node.states.clear();
		node.states.insert(states.begin(), states.end());// ASAP fix copying vector.
	}
}

bool DD::buildNextLayer(vector<int> &currentLayer, vector<int> &nextLayer, int index) {
	/*
	 * builds next layer from the given current layer.
	 * adds new child nodes and outgoing arcs to their respective maps.
	 * updates currentlayer's nodes and their arcs in the map.
	 * This function doesn't do reduction.
	 *
	 * return a boolean true if the next layer is exact, false otherwise.
	 */
	// TODO: reduction required?

	bool isExact = true;

	if (type == RESTRICTED) {

#if PRUNE == TRAIL
		{
			uint count = 0;

			for (const auto id: currentLayer) {

				DDNode &parentNode = nodes[id];
				const auto parentStates = parentNode.states;
				// INFO; states should contain -1.
				for (const auto decision: parentNode.states) {
					if (count >= MAX_WIDTH) {isExact = false; break;}

					DDNode node{lastInserted};
					DDArc arc{lastInserted, parentNode.id, node.id, decision};
					node.states = parentStates;
					if (decision != -1) node.states.erase(decision);
					//node.solutionVector = parentNode.solutionVector;
					//node.solutionVector.emplace_back(decision);
					node.incomingArcs.emplace_back(arc.id);
					parentNode.outgoingArcs.push_back(arc.id);
					// insert node and arc to map
					nodes.insert(std::make_pair(node.id, node));
					arcs.insert(std::make_pair(arc.id, arc));

					count++;
					lastInserted++;

					nextLayer.emplace_back(node.id);
				}
			}
		}

#elif PRUNE == RANDOM
	// remove nodes at random
#endif

	}
	else if(type == RELAXED){
#if PRUNE == TRAIL
		{
			int count = 0;
			int lastNodeId = 0;

			for (const auto id: currentLayer) {
				DDNode &parentNode = nodes[id];

				auto parentStates = parentNode.states;
				for (auto decision: parentNode.states) {

					if (count < MAX_WIDTH) {
						DDNode node{lastInserted};
						DDArc arc{lastInserted, id, node.id, decision};
						node.states = parentStates;
						if (decision != -1) node.states.erase(decision);
						//node.solutionVector = parentNode.solutionVector;
						//node.solutionVector.emplace_back(decision);
						node.incomingArcs.emplace_back(arc.id);
						parentNode.outgoingArcs.push_back(arc.id);
						// insert node and arc to map
						nodes.insert(std::make_pair(node.id, node));
						arcs.insert(std::make_pair(arc.id, arc));
						lastNodeId = node.id;
						nextLayer.emplace_back(node.id);
					} else { // create only arc and make it point to the last node.
						// QUESTION: do we create all arcs in relaxed version.
						DDArc arc{lastInserted, id, lastNodeId, decision};
						DDNode &node = nodes[lastNodeId];
						node.incomingArcs.emplace_back(arc.id);
						parentNode.outgoingArcs.emplace_back(arc.id);
						arcs.insert(std::make_pair(arc.id, arc));
						isExact = false;
					}
					count++;
					lastInserted++;
				}
			}
		}
#elif PRUNE == RANDOM
// merge nodes at random
#endif

	}
	return isExact;
}

//void ibuildNextLayer(vector<uint>& currentLayer, vector<uint>& nextLayer, int index) {
//	/*
//	 * build next layer from the given current layer.
//	 * should also create arcs and store them in the map.
//	 * should also create nodes and store them in the map.
//	 * updates both current layer and next layer
//	 */
//
//	//set<pair<set<int>,int>> allStates;
//
//	for (auto id: currentLayer){ // iterate through current layer.
//
//		auto& node = nodes[id];
//
//		// add zero state
//		DDNode zeroNode;
//		zeroNode.id = lastNode++; // TODO: need to add some more fields.
//		zeroNode.states = node.states;
//		DDArc zeroArc(lastArc++, zeroNode.id, id, 0);
//		node.outgoingArcs.push_back(zeroArc.id);
//		zeroNode.incomingArcs.push_back(zeroArc.id);
//		// insert to map.
//		nodes.insert({zeroNode.id, zeroNode});
//		arcs.insert({zeroArc.id, zeroArc});
//
//		nextLayer.push_back(zeroNode.id);
//
//		for (auto state: node.states){
//			// create arc and node
//			DDNode node1;// create a new node with some Id
//			node1.id = lastNode++;
//			auto states = node.states; // remove the state from states.
//			states.erase(std::remove(states.begin(), states.end(), state), states.end());
//			node1.states = std::move(states);
//			DDArc arc(lastArc++, node1.id, id, state);
//			node1.incomingArcs.push_back(arc.id);
//			node.outgoingArcs.push_back(arc.id); // NOTE: SHOULD BE ID OF THE CHILDARC, NOT CHILD NODE.
//			//DDArc arc(lastArc++, 10, 10, 10);
//			nodes.insert({node1.id, node1});
//			arcs.insert({arc.id, arc});
//			nextLayer.push_back(node1.id); // ASAP: remove the state from the list of states and
//		}
//	}
//
//	return std::move(nextLayer);
//}

void DD::mergeNodes(DDNode& node1, DDNode& node2) {
	/*
	 * Merge node2 with node1. Updates node1 attributes and removes node2 from dictionary.
	 * Used in reduceNodes().
	 */
	for (auto parent: node2.incomingArcs){ // update head id of the incoming arcs to node1.id
		auto& arc = arcs[parent];
		arc.head = node1.id; // TODO: change other attributes if needed
	}

	for (auto child: node2.outgoingArcs){ // update tail id of outgoing arcs to node1.id
		auto& arc = arcs[child];
		arc.tail = node1.id; // TODO: change other attributes if needed.
	}
	// update node1 incomingArcs to add node2's incomingArcs.
	node1.incomingArcs.insert(node1.incomingArcs.end(), node2.incomingArcs.begin(), node2.incomingArcs.end());
	node1.outgoingArcs.insert(node1.outgoingArcs.end(), node2.outgoingArcs.begin(), node2.outgoingArcs.end());

	node2.incomingArcs.clear();
	node2.outgoingArcs.clear();

}

void DD::reduceLayer(vector<int> &currentLayer) {
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

void DD::deleteArcById(int id){
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

	// delete node if no incoming arcs for the head node. TODO; INFO confirm with erfan later.
}

void DD::deleteNodeById(int id) {
	/*
	 * deletes node by its id.
	 * delete outgoing arcs vector. destructor deletes remaining attributes.
	 */
	auto& node = nodes[id];
	auto outgoingArcs = node.outgoingArcs;
	for (const auto arcId : outgoingArcs){
		deleteArcById(arcId);
	}
	// INFO not sure if this function shoudl remove the node from the  map.
	nodes.erase(id); // remove this if node deletion happens outside of the function.
}

inline DDNode DD::duplicate(const DDNode& node){
	// clone the node. for every outgoing arc, create new arc and point it to child node.
	DDNode dupNode(lastInserted++);
	dupNode.state2 = node.state2;
	dupNode.states = node.states;
	// ASAP how solution vector is copied
	// incoming arc is updated by the caller.

	for (const auto& outArcId: node.outgoingArcs){
		const auto& childNodeId = arcs[outArcId].head;
		DDArc newArc{lastInserted, dupNode.id, childNodeId, 0};
		dupNode.outgoingArcs.push_back(newArc.id);
		nodes[childNodeId].incomingArcs.push_back(newArc.id);
		arcs.insert(std::make_pair(newArc.id, newArc));
		lastInserted++;
	}

	return dupNode;
}

void DD::duplicateNode(int id){
	/*
	 * for every incoming arc, create new node (copy state parameters)
	 *  TODO: complete it today.
	 */
	auto& node = nodes[id];

	auto incomingArcs = node.incomingArcs;
	// skip first arc.
	for (uint i = 1; i < incomingArcs.size(); i++) {
		//
		int arcId = incomingArcs[i];
		DDArc& arc = arcs[arcId];

		// TODO compute state.
		bool feasible = true; // TODO: based on the cut, build if node is feasible.
		if (feasible){
			DDNode newNode = duplicate(node); // set incoming arc.
			// update incoming arc and make the arc to point to new node.
			newNode.incomingArcs.push_back(arcId);
			nodes.insert(std::make_pair(newNode.id, newNode));
			arc.head = newNode.id;
		}
		else {
			// delete arc.
			deleteArcById(arcId);
			auto& tailOutArcs = nodes[arc.tail].outgoingArcs;
			tailOutArcs.erase(std::find(tailOutArcs.begin(), tailOutArcs.end(), arcId));
			// remove arc
			arcs.erase(arcId);
		}
	}
	// remove the incoming arcs of original node.
	// compute state and keep it if node is feasible.
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

vi DD::solution() {
	/*
	 * Start from terminal node and iterate through incoming arcs and select arc.
	 *
	 */

	// change logic.

	// find maxVal parent node.
	int maxVal = 0;
	int maxNodeId = 0;
	vi path(tree.size()-1);
	int terminalNodeId = tree[tree.size()-1][0];
	const auto& terminalNode = nodes[terminalNodeId];

	for (const auto arcId: terminalNode.incomingArcs){ // find maxValue path.
		const auto& arc = arcs[arcId];
		if (arc.decision > maxVal){
			maxVal = arc.decision;
			maxNodeId = arc.tail;
		}
	}

	// maxVal and maxNodeId is populated.
	for (size_t l = tree.size()-2; l > 0; l--){ // check indexes later.
		const auto& node = nodes[maxNodeId];
		const auto& arc = arcs[node.incomingArcs[0]];
		path[l] = arc.weight;
		maxNodeId = arc.tail;
	}
	// prepend root solution to path.
	const auto& rootNode = nodes[tree[0][0]];
	vi finalPath(rootNode.solutionVector);
	finalPath.insert(finalPath.end(), path.begin(), path.end());
	return finalPath;

//
//	int maxVal = 0;
//	int maxNodeId = tree[tree.size()-1][0];
//	maxVal = nodes[maxNodeId].objVal;
//	vi path(tree.size());
//
//	for (size_t l = tree.size()-1; l > 0; l--){
//
//		for (const auto incomingId : nodes[maxNodeId].incomingArcs){
//			const auto& arc = arcs[incomingId];
//			const auto& parentNode = nodes[arc.tail];
//			if ( maxVal - arc.weight * 10 /* some coefficient */ == parentNode.objVal) {
//				// this node is parent.
//				maxNodeId = parentNode.id;
//				maxVal = maxVal - arc.weight; // ASAP multiply with coefficient.
//				path[l] = arc.weight;
//				break;
//			}
//		}
//	}
//	// copy the root's solution at front.
//	const auto& rootSolution = nodes[tree[0][0]].solutionVector;
//	for (int i = 0; i < rootSolution.size(); i++){
//		path[i] = rootSolution[i];
//	}

	return path;
}


vi DD::getSolutionVector(int nodeId){
	/*
	 * retrieve solution vector for a given node.
	 */

	auto& node = nodes[nodeId];
	int total = startTree + exactLayer-1;
	vi solutionVector(total);

	for (int i = exactLayer-1; i > 0; i--){
		// since all nodes upto exact layer has only parent
		const auto& inArc = arcs[node.incomingArcs[0]];
		solutionVector[startTree+i-1] =inArc.decision;
		node = nodes[inArc.tail];
	}
	// now node points to root node, add root's solution vector to it.
	for (int i = 0; i < node.solutionVector.size(); i++) solutionVector[i] = node.solutionVector[i];
	return solutionVector;
}

vector<DDNode> DD::getExactCutset() {
	/*
	 * returns a vector containing the nodes in the exact layer as cutset.
	 */
	//vector<DDNode> cutsetNodes(tree[exactLayer-1].size());
	vector<DDNode> cutsetNodes;
	int i = 0;
	for (const int id: tree[exactLayer-1]){
		auto node = nodes[id];
		node.solutionVector = getSolutionVector(id);
		//cutsetNodes[i] = node;
		cutsetNodes.push_back(node);
	}
	return cutsetNodes;
}
//int main(){
	// some processing order
	//vector<pair<int,int>> processingOrder = {{1,2},{2,3}};

	//DD diagram;

	//DDNode rootNode;
	//rootNode.states = unordered_set<int>({2,3,4});

	//diagram.build(processingOrder, rootNode, 0);

	//printTree(diagram.tree);
//}
