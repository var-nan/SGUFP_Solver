//
// Created by nandgate on 6/3/24.
//


#include "DD.h"

#include <chrono>
#include <execution>
#include <queue>
#include <limits>

/**
 * Compiles the decision diagram tree with given node as the root.
 */
optional<vector<Node_t>> DD::build(DDNode &node) {
	// node parameter should be initialized before calling this function. node should contain states.
	// the given node parameter will be inserted as the root of the tree.


	const auto& arcOrder = networkPtr->processingOrder;
	const auto& stateUpdateMap = networkPtr->stateUpdateMap;

	// set the node to the root.
	vector<ulint> currentLayer; // should be root layer.

	startTree = node.globalLayer; // LATER, size of solution vector of the root might be appropriate

	node.incomingArcs.clear();
	node.outgoingArcs.clear();
	node.id = 0; // make new node if necessary.
	node.nodeLayer = 0; // TODO change during branch and bound.
	nodes.insert(std::make_pair(node.id, node));
	currentLayer.push_back(node.id);
	// insert root layer to tree.
	tree.push_back(currentLayer);

	auto start = arcOrder.begin() + startTree;
	auto end = arcOrder.end();
	uint index = 0; // for node layer.

	for (; start != end; ++start){

		auto[a,b] = *start;
		if (stateUpdateMap.count(a)) { // update state of each node in the layer.
			const auto& newStates = stateUpdateMap.at(a);
			updateState(currentLayer, newStates);
		}
		//const unordered_set<int> temp = stateUpdateMap.at(a);
		vector<ulint> nextLayer;
		//nextLayer.reserve(MAX_WIDTH);
		isExact = buildNextLayer(currentLayer, nextLayer, networkPtr->hasStateChanged[a+1]);
		if (isExact) exactLayer++; // at last, this number should be exact layer number.
		//reduceLayer(nextLayer); // INFO not doing reduction.
		++index;
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}
	// terminal node layer.
	vector<ulint> terminalLayer;
	DDNode terminalNode {number.getNext()};
	terminalNode.nodeLayer = ++index;
	// current layer points to last layer of tree.
	for (const auto& id: currentLayer){
		// only add single arc for each node in the last layer.
		ulint arcId = number.getNext();
		auto& parentNode = nodes[id];
		// create arc add to incoming of terminal node.
		DDArc arc{arcId, id, terminalNode.id, 1};
		arc.weight = std::numeric_limits<double>::max();
		terminalNode.incomingArcs.push_back(arcId);
		parentNode.outgoingArcs.push_back(arcId);
		arcs.insert(make_pair(arcId, arc));
	}

	nodes.insert(make_pair(terminalNode.id, terminalNode));
	terminalLayer.push_back(terminalNode.id);
	tree.push_back(terminalLayer);

	// generate cutset.
	//cutSet = generateExactCutSet(); TODO uncomment this.


	#ifdef DEBUG
	// displayArcLabels();
		displayStats();
	#endif
	if (exactLayer == index-1) return {};
	return generateExactCutSet();
	// return (exactLayer == index+1) ? generateExactCutSet() : std::nullopt;
}

inline void DD::updateState(const vector<ulint> &currentLayer, const set<int> &states){

	for (auto id: currentLayer){
		auto& node = nodes[id];
		node.states.clear();
		node.states.insert(states.begin(), states.end());
	}
}

/**
 * Builds the next layer given the current layer corresponding to the type of DD.
 * Strategies to build the next layer can be modified by changing '_STRATEGY' macro.
 * Returns boolean true if the newly built layer is an exact layer, else returns boolean false.
 */
bool DD::buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, bool stateChangesNext) {
	/*
	 * builds next layer from the given current layer.
	 * adds new child nodes and outgoing arcs to their respective maps.
	 * updates current layer's nodes and their arcs in the map.
	 * This function uses reduction only during compilation of relaxed tree.
	 */

	bool isExact = true;

	if (type == RESTRICTED) {

		#if RESTRICTED_STRATEGY == 1
		{
			if (currentLayer.size() >= MAX_WIDTH) { // new strategy after tree becomes non-exact.
				buildNextLayer6(currentLayer, nextLayer);
				return false;
			}
			uint count = 0;

			for (const auto id: currentLayer) {

				DDNode &parentNode = nodes[id];
				const auto parentStates = parentNode.states;
				// INFO; states should contain -1.
				for (auto start = parentNode.states.rbegin(); start != parentNode.states.rend(); ++start){
					auto decision = *start; // iterate reverse order.
				// for (const auto decision: parentNode.states) {
					if (count >= MAX_WIDTH) {isExact = false; break; } // jump with goto, instead of this.
					auto lastInserted = number.getNext();
					DDNode node{lastInserted};
					DDArc arc{lastInserted, parentNode.id, node.id, decision};
					node.states = parentStates;
					if (decision != -1) node.states.erase(decision);
					//node.solutionVector = parentNode.solutionVector; // solutions are computed during bulding cutset.
					//node.solutionVector.emplace_back(decision);
					node.incomingArcs.emplace_back(arc.id);
					node.nodeLayer = parentNode.nodeLayer+1;
					node.globalLayer = parentNode.globalLayer+1;
					parentNode.outgoingArcs.push_back(arc.id);
					// insert node and arc to map
					nodes.insert(std::make_pair(node.id, node));
					arcs.insert(std::make_pair(arc.id, arc));
					count++;

					nextLayer.emplace_back(node.id);
				}
			}
		}

		#elif RESTRICTED_STRATEGY == 2
			// remove nodes at random
		#endif

	}
	else if(type == RELAXED){
		/*
		 * RELAXED_STRATEGY . make relaxed arcs to point to last inserted node in
		 * the current layer.
		 */
		#if RELAXED_STRATEGY == 3
		{
			int count = 0;
			ulint lastNodeId = 0;

			for (const auto id: currentLayer) {
				DDNode &parentNode = nodes[id];

				auto parentStates = parentNode.states;
				for (auto decision: parentNode.states) {

					auto lastInserted = number.getNext();
					if (count < MAX_WIDTH) {
						DDNode node{lastInserted};
						DDArc arc{lastInserted, id, node.id, decision};
						node.states = parentStates;
						if (decision != -1) node.states.erase(decision);
						//node.solutionVector = parentNode.solutionVector;
						//node.solutionVector.emplace_back(decision);
						node.incomingArcs.emplace_back(arc.id);
						parentNode.outgoingArcs.push_back(arc.id);
						node.nodeLayer = index;
						// insert node and arc to map
						nodes.insert(std::make_pair(node.id, node));
						arcs.insert(std::make_pair(arc.id, arc));
						lastNodeId = node.id;
						nextLayer.emplace_back(node.id);
					} else { // create only arc and make it point to the last node.
						// QUESTION: do we create all arcs in relaxed version.
						DDArc arc{lastInserted, id, lastNodeId, decision};
						DDNode &node = nodes[lastNodeId];
						node.states.insert(parentStates.begin(), parentStates.end()); // insert parent states
						if(decision != -1) node.states.erase(decision); // remove decision
						node.incomingArcs.emplace_back(arc.id);
						parentNode.outgoingArcs.emplace_back(arc.id);
						arcs.insert(std::make_pair(arc.id, arc));
						isExact = false;
					}
					count++;
					//lastInserted++;
				}
			}
		}
		//#endif

		#elif RELAXED_STRATEGY == 2
		{
			// merge all outgoing arcs of parent to the same children.

			int count = 0;
			ulint lastNode;

			if (currentLayer.size() >= MAX_WIDTH) {
				// for each parent, create only one child.
				for (const auto& id: currentLayer) {
					// create child node.
					auto& parentNode = nodes[id];
					auto parentStates = parentNode.states;

					auto lastInserted = number.getNext();
					DDNode node{lastInserted};
					node.states = parentStates;

					for (auto decision : parentStates) {
						//lastInserted = number.getNext();
						DDArc arc{lastInserted, id, node.id, decision};
						//node.states = parentStates;
						//if (decision != -1) node.states.erase(decision); // INFO not required, because relaxation.
						parentNode.outgoingArcs.push_back(lastInserted);
						node.incomingArcs.push_back(lastInserted);

						arcs.insert(std::make_pair(lastInserted, arc));
						lastInserted = number.getNext();
						//nodes.insert(std::make_pair(lastInserted, node));
					}
					// insert node to the map and node id to tree.
					nodes.insert(std::make_pair(node.id, node));
					nextLayer.push_back(node.id);
					isExact = false; // TODO is this correct?
				}
			}
			else { // build complete layer. this might cross MAX_WIDTH.
				for (const auto& id : currentLayer) {
					auto& parentNode = nodes[id];
					auto parentStates = parentNode.states;

					for(auto decision : parentStates) {
						auto nextId = number.getNext();
						DDNode node{nextId};
						DDArc arc{nextId, id, node.id, decision};
						node.states = parentStates;
						if (decision != -1) node.states.erase(decision);
						node.incomingArcs.push_back(arc.id);
						parentNode.outgoingArcs.push_back(arc.id);

						nodes.insert(std::make_pair(node.id, node));
						arcs.insert(std::make_pair(arc.id, arc));

						nextLayer.push_back(node.id);
					}
				}
				if (nextLayer.size() > MAX_WIDTH) isExact = false;
			}
		}
		#elif RELAXED_STRATEGY == 1

		#endif
	}
	else { // exact tree.
		#if EXACT_STRATEGY == 1
		{
			 // build complete tree with state reduction.
			 if (stateChangesNext) { // next DD layer will undergo state change, so create only one node.

			 	// if current layer is trouble maker, then handle edge case.
			 	const auto& layer = nodes[currentLayer.front()].globalLayer;
			 	bool res = stateChangesNext; //networkPtr->troubleMaker[layer] == 1;

				 auto nextNodeId = number.getNext();
				 DDNode newNode{nextNodeId};
				 //newNode.nodeLayer = index;

				 for (const auto id : currentLayer) {
					 assert(nodes.count(id));
					 auto &node = nodes[id];
				 	 newNode.nodeLayer = node.nodeLayer+1;
				 	 newNode.globalLayer = node.globalLayer+1;
					 for (auto decision: node.states) {
					 	if (res && node.states.size() > 1 && decision == -1) {
					 		// do not create arc.
					 		// cout << "Trouble maker arc is not created." << endl;
					 		continue;
					 	}
						 auto nextId = number.getNext();
						 DDArc newArc{nextId, id, nextNodeId, decision};
						 node.outgoingArcs.push_back(nextId);
						 newNode.incomingArcs.push_back(nextId);
						 arcs.insert(make_pair(nextId, newArc));
					 }
				 } // state of this node is updated in the build() in next iteration.
				 nextLayer.push_back(nextNodeId);
				 nodes.insert(make_pair(nextNodeId, newNode));
			 }
			 else {
				 vector<DDNode> nodesVector;
				 unordered_set<tuple<set<int>, int, int>, tuple_hash, tuple_equal> allStates;
				 int j = 0;

				 for (const auto id: currentLayer) {
					 assert(nodes.count(id));
					 auto &node = nodes[id];
					 auto statesCopy = node.states;

					 for (auto decision: node.states) {
						 auto newStates(statesCopy);
						 if (decision != -1) newStates.erase(decision);
						 auto nextId = number.getNext();

						 // if newStates already in allStates, update exising node in nodesVector.
						 // auto [it, isInserted] = allStates.insert({newStates, 0, j});
						 // if (isInserted) { // create new Node
							 DDNode newNode{nextId};
							 DDArc newArc{nextId, id, nextId, decision};
							 newNode.nodeLayer = node.nodeLayer+1;
						 	 newNode.globalLayer = node.globalLayer+1;
							 newNode.states = newStates;
							 newNode.incomingArcs.push_back(nextId);
							 node.outgoingArcs.push_back(nextId);
							 arcs.insert(make_pair(nextId, newArc));
							 nodesVector.push_back(newNode);
							 j++;
						 // } else { // state already exists, update existing node in nodesVector.
							//  auto [tempState, state2, pos] = *(it);
							//  auto &prevNode = nodesVector[pos];
							//  DDArc newArc{nextId, id, prevNode.id, decision};
							//  prevNode.incomingArcs.push_back(nextId);
							//  node.outgoingArcs.push_back(nextId);
							//  arcs.insert(make_pair(nextId, newArc));
						 // }
					 }
				 }

				 // populate nextLayer. LATER. move nodes from vector.
				 for (auto& node : nodesVector) {
					 nextLayer.push_back(node.id);
					 nodes.insert(make_pair(node.id, node));
				 }
			 }
		}
		#elif EXACT_STRATEGY == 2
			// without state reduction.
		#endif
	}
	return isExact;
}

/**
 * Creates a child node for the given parent node, and the corresponding arc
 * between them. Inserts the child node and the arc to the map. Returns the
 * id of the child node, this id is same for the arc between them.
 * @param parent parent node
 * @param decision decision of the arc between parent and child node.
 * @return id of the child node (= id of arc between parent and child).
 */
ulint DD::createChild(DDNode& parent, int decision) {
	// create child for this node with the decision.
	auto nextId = number.getNext();
	DDNode newNode{nextId};
	newNode.states = parent.states;
	if (decision != -1) newNode.states.erase(decision);
	DDArc arc {nextId, parent.id, nextId, decision};
	newNode.incomingArcs.push_back(nextId);
	parent.outgoingArcs.push_back(nextId);
	nodes.insert(make_pair(nextId, newNode));
	arcs.insert(make_pair(nextId, arc));
	return nextId;
}

/**
 * Build next layer with minimum number of nodes that can be created with -1.
 * Assumes next layer is not exact. DO NOT CALL this function if DD is still exact.
 * @param currentLayer
 * @param nextLayer
 */
void DD::buildNextLayer2(vector<ulint> &currentLayer, vector<ulint> &nextLayer) {
	uint count = 0;

	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		if (node.states.size() > 1) {
			for (auto decision: node.states) {
				if (decision == -1) continue;
				auto childId = createChild(node, decision);
				nextLayer.push_back(childId);
				if (count++ >= MAX_WIDTH) return;
			}
		}
		else {
			auto childId = createChild(node, -1);
			nextLayer.push_back(childId);
			if (count++ >= MAX_WIDTH) return;
		}
	}
}


void DD::buildNextLayer3(vector<ulint>& currentLayer, vector<ulint>& nextLayer) {
	uint count = 0;

	// for each node, create only one child
	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		if (node.states.size() > 1) {
			for (auto state: node.states) {
				if (state == -1) continue;
				// create child
				auto childId = createChild(node, state);
				nextLayer.push_back(childId);
				break;
			}
		}
		else {
			auto childId = createChild(node, -1);
			nextLayer.push_back(childId);
		}
		count++;
		if (count >= MAX_WIDTH) return;
	}
}

void DD::buildNextLayer4(vector<ulint> &currentLayer, vector<ulint> &nextLayer) {
	// select arc with maximum reward. for each parent create single child

	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		int decision = -1;
		if (node.states.size() > 1) // select arc with maximum reward.
			decision = networkPtr->getBestArc(node.states);
		auto childId = createChild(node, decision);
		nextLayer.push_back(childId);
	}
}

void DD::buildNextLayer5(vector<ulint> &currentLayer, vector<ulint> &nextLayer) {
	// match out arc that has max reward to the current incoming arc.

	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		// get layer number
		auto inArcReward = networkPtr->layerRewards[node.globalLayer];
		// get best arc
		auto bestState = networkPtr->getBestArc(node.states);
		int decision = -1;
		if (bestState != -1 && networkPtr->networkArcs[bestState].rewards[0] + inArcReward >= 0) {
			// match this out arc with in arc.
			decision = bestState;
		}
		auto childId = createChild(node, decision);
		nextLayer.push_back(childId);
	}
}

void DD::buildNextLayer6(vector<ulint> &currentLayer, vector<ulint> &nextLayer) {
	// if current layer is trouble maker, then handle troublemaker node.
	const auto& layerNumber = nodes[currentLayer.front()].globalLayer;
	// bool res = (networkPtr->troubleMaker[layerNumber] == 1);
	// if (networkPtr->troubleMaker[layerNumber] == 1) {
	if (networkPtr->hasStateChanged[layerNumber+1]){ // state changes in next layer, thus remove -1 arcs from this layer.
		// for all the nodes, remove paths that contains all -1's for this node.
		for (const auto id: currentLayer) {
			auto& node = nodes[id];
			int decision = -1;
			for (auto start = node.states.rbegin(); start != node.states.rend(); ++start) {
				auto state = *start;
				if (state > decision) decision = state;
			}
			// for (auto state: node.states){ if (state > decision) decision = state;}
			auto childId = createChild(node, decision);
			nextLayer.push_back(childId);
		}
		return;
	}
	// select one state from states at random
	srand(time(nullptr));
	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		vi states{node.states.begin(), node.states.end()};
		auto s = rand()%states.size();
		auto decision = states.back();
		auto childId = createChild(node, decision);
		nextLayer.push_back(childId);

	}
}

/**
 * Merges the second DDNode into the first DDNode.
 * Updates fist node's attributes and deletes second node from the map.
 */
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

void DD::reduceLayer(vector<ulint> &currentLayer) {
	/*
	 * Merge two nodes that containing same state. Update incoming and outgoing arcs to point to merged node.
	 */

	queue<ulint> q;
	unordered_set<tuple<set<int>,int,int>,tuple_hash,tuple_equal> allStates;
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
			q.push(i);
		}
		j++;
	}
	queue tempQ {q};
	while (!tempQ.empty()) {number.setNext(tempQ.front()); tempQ.pop();}
	// remove filtered nodeIds from current layer.
	currentLayer.erase(std::remove_if(currentLayer.begin(), currentLayer.end(),
	                                  [&q](const int x) mutable{
												if ((!q.empty()) && (q.front() == x)) { q.pop(); return true;}
												return false;
											}),currentLayer.end());

}

/**
 * Deletes the arc from the arcs map given its id.
 */
void DD::deleteArcById(ulint id){
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

/**
 *
 */
void DD::deleteNodeById(ulint id) {
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
	deletedNodeIds.insert(id); // this id should be removed from the tree layer, after refinement.
}

inline DDNode DD::duplicate(const DDNode& node){
	// clone the node. for every outgoing arc, create new arc and point it to child node.
	auto lastInserted = number.getNext();
	DDNode dupNode(lastInserted);
	dupNode.state2 = node.state2;
	dupNode.states = node.states;
	// LATER not sure how to copy solution vector
	// incoming arc is updated by the caller.

	for (const auto& outArcId: node.outgoingArcs){
		const auto& childNodeId = arcs[outArcId].head;
		lastInserted = number.getNext();
		DDArc newArc{lastInserted, dupNode.id, childNodeId, arcs[outArcId].decision};
		dupNode.outgoingArcs.push_back(newArc.id);
		nodes[childNodeId].incomingArcs.push_back(newArc.id);
		arcs.insert(std::make_pair(newArc.id, newArc));
		//lastInserted++;
	}

	return dupNode;
}

void DD::duplicateNode(ulint id){
	/*
	 * for every incoming arc, create new node (copy state parameters)
	 *  TODO: complete it today.
	 */
	auto& node = nodes[id];

	auto incomingArcs = node.incomingArcs;
	// skip first arc.
	for (uint i = 1; i < incomingArcs.size(); i++) {
		//
		ulint arcId = incomingArcs[i];
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


vi DD::solution() {
	/*
	 * Start from terminal node and iterate through incoming arcs and select arc.
	 * Returns the maximum path from the tree.
	 */

	// find maxVal parent node.
	double maxVal = std::numeric_limits<double>::lowest();
	ulint maxNodeId = 0;

	vi path;
	path.reserve(tree.size()-1);
	ulint terminalNodeId = tree[tree.size()-1][0];
	const auto& terminalNode = nodes[terminalNodeId];

	for (const auto arcId: terminalNode.incomingArcs){ // find maxValue path.
		const auto& arc = arcs[arcId];
		if (arc.weight > maxVal){
			maxVal = arc.weight; // change to objective value later.
			maxNodeId = arc.tail;
		}
	}

	// maxVal and maxNodeId is populated.
	for (size_t l = tree.size()-2; l > 0; l--){
//		const auto& node = nodes[maxNodeId];
//		const auto& arc = arcs[node.incomingArcs[0]];
//		path[l] = arc.decision;
//		maxNodeId = arc.tail;
		const auto& node = nodes[maxNodeId];
		for (const auto& arcId : node.incomingArcs){
			const auto& arc = arcs[arcId];
			const auto& tailNode = nodes[arc.tail];

			if ((tailNode.state2 + arc.weight) == node.state2){
				maxNodeId = tailNode.id;
				path.push_back(arc.decision);
				break;
			}
		}
	}
	// prepend root solution to path.
	const auto& rootNode = nodes[tree[0][0]];
	vi finalPath(rootNode.solutionVector);
	finalPath.insert(finalPath.end(), path.rbegin(), path.rend());
	return finalPath;
}

/**
 * Returns solution vector for a given node in the exact layer. function should only used for
 * a node in the exact cutset.
 * @param nodeId - id of node in the cutset.
 * @return - vector of solutions from root to the given node.
 */
vi DD::computePathForExactNode(ulint nodeId) const{
	DDNode node = nodes.at(nodeId);
	vi solutionVector;
	// iterate through all nodes from current node till root.
	while (!node.incomingArcs.empty()) {
		const auto& inArc = arcs.at(node.incomingArcs[0]);
		solutionVector.push_back(inArc.decision);
		node = nodes.at(inArc.tail);
	}
	// add root's solution.
	vi result(node.solutionVector);
	result.insert(result.end(), solutionVector.rbegin(), solutionVector.rend());
	return result;
}

/**
 * Creates a copy of given node. copies only the id,states
 * and node-layer attributes.
 * @param node
 * @return
 */
DDNode copyNode(const DDNode& node) {
	DDNode newNode{node.id};
	newNode.states = node.states;
	newNode.nodeLayer = node.nodeLayer;
	return newNode;
}

/**
 * Computes and returns a vector of nodes that are in exact cutset.
 * @return vector of nodes that are in exact cutset.
 */
vector<Node_t> DD::generateExactCutSet() const {
	/*
	 * returns a vector containing the nodes in the exact layer as cutset.
	 */
	vector<Node_t> cutsetNodes;
	cutsetNodes.reserve(MAX_WIDTH);

	for (const auto id: tree[exactLayer]){
		const auto& node = nodes.at(id);
		vi states{node.states.begin(), node.states.end()};
		cutsetNodes.emplace_back(states, computePathForExactNode(id),
						std::numeric_limits<double>::lowest(),std::numeric_limits<double>::max(),
						node.globalLayer);
	}
	return cutsetNodes;
}

vector<Node_t> DD::getExactCutSet() {
	return cutSet;
}


static vector<double> helperFunction(const Network &network, const Cut &cut) {
	/*
	 * generates the lower bounds for each DD layer. used in feasibilitycut refinement.
	 */
	vector<double> lowerBounds(network.processingOrder.size());
	const auto& coeff = cut.cutCoeff;

	for (const auto[id_t, arcId]: network.processingOrder){
		const auto& arc = network.networkArcs[arcId];
		auto i = arc.tailId;
		auto q = arc.headId;

		const auto& node = network.networkNodes[q];
		double max = 0;
		for (const auto j: node.outNodeIds){
			auto c = coeff.at(make_tuple(static_cast<int>(i),static_cast<int>(q),static_cast<int>(j)));
			if (c > max) max = c;
		}
		lowerBounds[id_t] = max;
	}
	// compute suffix sum
	double pref = 0;
	for (auto start = lowerBounds.rbegin(); start != lowerBounds.rend(); start++){
		*start = *start + pref;
		pref = *start;
	}
	return lowerBounds; // automatic copy elision?
}

void DD::refineTree(const Network &network, Cut cut) {

	if (cut.cutType == OPTIMALITY){
		applyOptimalityCut(network, cut);
	}
	else applyFeasibilityCut(network, cut);
}


void DD::applyFeasibilityCut(const Network &network, const Cut &cut) {

	vector<double> lowerBounds = helperFunction(network, cut);

	const auto& coeff = cut.cutCoeff;
	double RHS = cut.RHS;

}

void DD::applyOptimalityCut(const Network &network, const Cut &cut) {
	// TODO: incorporate semiroot partial solution.

}

void DD::deleteArc(DDNode& tailNode, DDArc& arc, DDNode& headNode){
	/*
	 * deletes the given arc from the arcs map.
	 * updates the head node's incomingArcs vector and tail node's outgoingArcs vector.
	 */
	assert(tailNode.id == arc.tail && headNode.id == arc.head);
	assert(arcs.count(arc.id));
	auto pos = std::find(headNode.incomingArcs.begin(), headNode.incomingArcs.end(), arc.id);
	if (pos != headNode.incomingArcs.end()) headNode.incomingArcs.erase(pos);
	auto pos2 = std::find(tailNode.outgoingArcs.begin(), tailNode.outgoingArcs.end(), arc.id);
	if (pos2 != tailNode.outgoingArcs.end()) tailNode.outgoingArcs.erase(pos2);
	arcs.erase(arc.id);
}

/**
 * Deletes the given node from the nodes container. Inserts the node id to the list of
 * deleted node Ids.
 *
 * A call to this function should be made iff both incoming and outgoing arcs are empty.
 */
void DD::deleteNode(DDNode& node){
	assert(nodes.count(node.id));
	deletedNodeIds.insert(node.id);
	nodes.erase(node.id);
}

/**
 * Deletes the subtree of the given node recursively.
 *
 * Removes only the child nodes that have one incoming arc.
 */
void DD::topDownDelete(ulint id) { // hard delete function.
	/*
	 * start from the current node and recursively delete the children nodes
	 * until next node is terminal node or node with multiple incoming arcs.
	 * remove the last outgoing arc of the last node.
	 */
	{
		assert(nodes.count(id));
		auto &node = nodes[id];
		// remove the incoming arc here.
		//deleteArcById(node.incomingArcs[0]);
		//assert(!node.outgoingArcs.empty());
		auto arcsToDelete = node.outgoingArcs;

		for (auto outArcId: arcsToDelete) {
			// for each outArcId, find its head and apply topDownDelete() if head has single incoming arc.
			assert(arcs.count(outArcId));
			auto& outArc = arcs.at(outArcId);
			auto childId = outArc.head;
			assert(nodes.count(childId));
			auto &childNode = nodes[childId];
			deleteArc(node, outArc, childNode);
			if (childNode.incomingArcs.empty()) topDownDelete(childId); // orphan node
		}
		assert(node.outgoingArcs.empty());
		deleteNode(node);
	}
	// here, node has no outgoing arcs and incoming arcs left,
	// also remove the node id from the tree layer.
	//deleteNode(node);
}

/**
 * Removes the node and its associated nodes from the nodes map and updates the
 * tree vector if deletion is not occurring in batch.
 * @param id Id of the node to be removed.
 * @param isBatch updates the tree vector if parameter is false. Default value is false.
 */
void DD::removeNode(ulint id, bool isBatch){
	/*
	 * NOTE:
	 */
	// if current node has multiple parents and multiple children, do this.
	auto& node = nodes[id];
	auto outArcs = node.outgoingArcs;
	assert(!node.outgoingArcs.empty());
	for (auto childArcId : outArcs){
		assert(arcs.count(childArcId));
		auto& childArc = arcs[childArcId];
		assert(nodes.count(childArc.head));
		auto& child = nodes[childArc.head];
		deleteArc(node, childArc, child);
		if (child.incomingArcs.empty()) topDownDelete(child.id); // orphan node
	}
	//if (node.incomingArcs.size() > 1) { // multiple parents
	auto incomingArcs = node.incomingArcs;
	assert(!node.incomingArcs.empty());
	for (auto arcId: incomingArcs){
		assert(arcs.count(arcId));
		auto& arc = arcs[arcId];
		assert(nodes.count(arc.tail));
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		if (parentNode.outgoingArcs.empty()) bottomUpDelete(parentNode.id); // orphan node
	}
	assert(node.incomingArcs.empty() && node.outgoingArcs.empty());
	// actual delete.
	deleteNode(node);
	assert(!nodes.count(id));

	if (!isBatch) updateTree();	// no batch deletion.
}

/**
 * Removes the deleted node ids from the tree. Resets the deletedNodeIds variable after
 * updating the tree. Adds the deletedNodeIds to the number.
 *
 * Optimized to call the function when removing nodes in batch.
 */
void DD::updateTree() {
	int n_removed = deletedNodeIds.size();
	auto& deletedNodeIds_l = this->deletedNodeIds;
	auto f = [&n_removed, &deletedNodeIds_l] (ulint x) mutable {
		if (deletedNodeIds_l.count(x)) { n_removed--; return true;}
		return false;
	};
	#ifdef DEBUG
		// cout << "Removing " << n_removed << " nodes from tree." << endl;
	#endif
	for (int i = tree.size()-2; i > 0; i--){ // iterate every layer until all deleted Ids are removed from tree.
		if (n_removed){
			auto& layer = tree[i];
			layer.erase(std::remove_if(layer.begin(), layer.end(), f), layer.end());
		}
		else break;
	}
	// add deleted Ids to the numbers.
	auto& num = number;
	std::for_each(deletedNodeIds.begin(), deletedNodeIds.end(), [&num](ulint x) mutable{num.setNext(x);});
	deletedNodeIds.clear(); // clear the deleted NodeIds set.
}

/**
 * Removes the given node ids from the DD and updates the tree at once.
 * @param ids vector of node ids to be removed.
 */
void DD::batchRemoveNodes(const vulint& ids) {
	// remove all the ids without updating the tree.
	for (auto id : ids) removeNode(id, true);
	// update tree
	updateTree();
}


void DD::bottomUpDelete(ulint id){
	// INFO this node might contain multiple incoming parents, but might contain one (or zero) children.

	// for each incoming arc, remove arc and call soft delete on incoming node (iff has single outgoing arc).
	auto& node = nodes[id];
	auto incomingArcs = node.incomingArcs;
	for (const auto& arcId : incomingArcs){
		auto& arc = arcs[arcId];
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		//deleteArcById(arcId); // delete this incoming arc.
		if (parentNode.outgoingArcs.empty()){
			bottomUpDelete(parentNode.id);
		}
	}
	// delete this node from nodes.
	deleteNode(node);

	// phase 1: recursively reach parents if they don't have multiple children.
	// if multiple incoming arcs? all bottomUpDelete () on each of the incoming nodes.
	// phase 2: recursively reach down the tree until reaching terminal node or children with multiple incoming arcs.
	// delete each node and arc in between.

}

void DD::applyFeasibilityCutRestricted(const Network &network, const Cut &cut) {

	// set the root state to RHS. // TODO: during refinement, compute new RHS for the subroot.

	nodes[0].state2 = cut.RHS;
	// heuristic to remove nodes.
	vector<double> lowerBounds = helperFunction(network, cut);
	lowerBounds.push_back(0);

	// set the state2 to all nodes to zero.

	for (size_t layerNo = 1; layerNo < tree.size()-1; layerNo++){

		auto netArc = network.processingOrder[layerNo-1].second;
		auto i_NetNodeId = network.networkArcs[netArc].tailId;
		auto q_NetNodeId = network.networkArcs[netArc].headId;

		vector<ulint> IdsToBeRemoved{};

		for(auto id: tree[layerNo]){
			auto& node = nodes[id];
			auto inArc = node.incomingArcs[0];
			auto& arc = arcs[inArc];
			auto parentId = arcs[inArc].tail;

			if (arc.decision != -1){
				auto j_NetNodeId = network.networkArcs[arc.decision].headId;
				arc.weight = cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId));
				auto newState = nodes[parentId].state2 + arc.weight;
				if (newState + lowerBounds[layerNo] >= 0) node.state2 = newState;
				else IdsToBeRemoved.push_back(id);
			}
			else {
				node.state2 = nodes[parentId].state2;
				if (node.state2 + lowerBounds[layerNo] < 0)
					IdsToBeRemoved.push_back(id);
			}
		}
		// remove the nodes
		batchRemoveNodes(IdsToBeRemoved);
	}
}

void DD::applyFeasibilityCutRelaxed(const Network &network, const Cut &cut) {
	// nodes might contain multiple parents.

	nodes[0].state2 = cut.RHS; // TODO: during branch and bound, compute new RHS for the subroot.
	// compute heuristics
	vector<double> lowerBounds = helperFunction(network, cut);
	lowerBounds.push_back(0);

	// terminal node
	auto& terminalNode = nodes[tree[tree.size()-1][0]];

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArc = network.processingOrder[layer-1].second; // TODO during branch and bound, compute correct index.
		auto i_NetNodeId = network.networkArcs[netArc].tailId;
		auto q_NetNodeId = network.networkArcs[netArc].headId;

		vector<ulint> nodesToRemove; // deleted nodes, remove ids from the layer.
		vector<ulint> nodesToAdd; // duplicated nodes.

		for (auto id: tree[layer]) {
			auto& node = nodes[id];
			node.state2 = 0;
			node.nodeLayer = layer; // TODO remove this?
			// for all arcs (except first), duplicate node if their state is feasible
			auto incomingArcs = node.incomingArcs;

			if (incomingArcs.size() > 1) {
				// duplicate nodes and outgoing arcs
				vulint arcsToRemove;
				for (size_t i = 1; i < incomingArcs.size(); i++) {
					auto inArc = incomingArcs[i];
					auto& arc = arcs[inArc];
					const auto& parentNode = nodes[arc.tail];

					auto j_NetNodeId = (arc.decision != -1) ? network.networkArcs[arc.decision].headId : 0;
					arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
					double newState = parentNode.state2 + arcs[inArc].weight;

					if (newState + lowerBounds[layer] >= 0){
						auto nextIndex = number.getNext();
						DDNode newNode{nextIndex};
						newNode.state2 = newState;
						newNode.nodeLayer = node.nodeLayer;
						arc.head = nextIndex;
						// if states changed in current layer, copy states from duplicating node, else copy from parent node.
						if (network.hasStateChanged[newNode.nodeLayer]) newNode.states = network.stateUpdateMap.at(layer);
						else { // copy states of parent, and remove the decision if decision is not -1
							newNode.states = parentNode.states;
							if (arc.decision != -1) newNode.states.erase(arc.decision);
						}
						newNode.incomingArcs = {arc.id};
						// handle pre-terminal layer separately.
						if (newNode.nodeLayer != tree.size()-2) {
							for (auto outArcId: node.outgoingArcs) {
								const auto &outArc = arcs[outArcId];
								if (!newNode.states.count(outArc.decision)) continue;
								auto nextArcId = number.getNext();
								DDArc newArc{nextArcId, nextIndex, outArc.head, outArc.decision};
								arcs.insert(make_pair(nextArcId, newArc));
								newNode.outgoingArcs.push_back(nextArcId);
								nodes[newArc.head].incomingArcs.push_back(nextArcId);
							}
						}
						else { // layer above terminal.
							const auto& outArc = arcs[node.outgoingArcs[0]];
							auto nextArcId = number.getNext();
							DDArc newArc{nextArcId, newNode.id, terminalNode.id, outArc.decision };
							newArc.weight = outArc.weight;
							terminalNode.incomingArcs.push_back(nextArcId);
							newNode.outgoingArcs.push_back(nextArcId);
							arcs.insert(make_pair(nextArcId, newArc));
						}

						nodes.insert(make_pair(nextIndex, newNode));
						nodesToAdd.push_back(nextIndex);
					}
					else arcsToRemove.push_back(arc.id);
				}
				for (const auto badArcId: arcsToRemove) deleteArcById(badArcId);
			}
			// process first arc. at this point, the node should have one incoming arc.
			auto inArc = incomingArcs[0];
			auto& arc = arcs[inArc];
			node.incomingArcs = {inArc};
			const auto& parentNode = nodes[arcs[inArc].tail];

			if (!network.hasStateChanged[layer]) {
				node.states = parentNode.states;
				if (arc.decision != -1) node.states.erase(arc.decision);
			}
			else node.states = network.stateUpdateMap.at(layer);

			auto j_NetNodeId = (arc.decision != -1) ? network.networkArcs[arc.decision].headId : 0;
			arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple( i_NetNodeId, q_NetNodeId, j_NetNodeId )) : 0;
			auto newState = parentNode.state2 + arc.weight;
			if (newState + lowerBounds[layer] >= 0) {
				node.state2 = newState;
				// remove arcs that are found to be infeasible (not present in the node's states).
				if (node.nodeLayer != tree.size()-2) {
					auto outArcs = node.outgoingArcs;
					for (auto outArc: outArcs) {
						if (!node.states.count(arcs[outArc].decision)) deleteArcById(outArc);
					}
				}
			}
			else nodesToRemove.push_back(id);

		}
		// remove nodes that are marked for deletion.
		batchRemoveNodes(nodesToRemove);
		// insert new Ids to current layer.
		tree[layer].insert(tree[layer].end(), nodesToAdd.begin(), nodesToAdd.end());
	}
}

/**
 * Applies given optimality cut to the tree starting with non-root layer until the
 * terminal layer. The optimality cut might prune some solutions that are infeasible
 * to some scenarios. Thus, the function may remove nodes that are infeasible based
 * on the cut.
 *
 * @param cut 
 */
void DD::applyOptimalityCutRestricted(const Cut &cut) {
	/*
	 * 
	 */

	nodes[0].state2 = cut.RHS; // TODO change this during branch and bound.

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArc = networkPtr->processingOrder[layer-1].second;
		auto i_NetNodeId = networkPtr->networkArcs[netArc].tailId;
		auto q_NetNodeId = networkPtr->networkArcs[netArc].headId;

		for (auto id: tree[layer]) {
			auto& node = nodes[id];
			auto inArcId = node.incomingArcs[0];
			auto& arc = arcs[inArcId];
			auto& parentNode = nodes[arc.tail];

			auto j_NetNodeId = (arc.decision != -1) ? networkPtr->networkArcs[arc.decision].headId : 0;
			arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
			node.state2 = arc.weight + parentNode.state2;
		}
	}

	// Update the state and arc weights of the terminal.
	auto& terminalNode = nodes[tree[tree.size()-1][0]];
	double terminalState = std::numeric_limits<double>::lowest();

	for (const auto arcId : terminalNode.incomingArcs){
		auto& arc = arcs[arcId];
		auto newState = nodes[arc.tail].state2;
		arc.weight = (newState < arc.weight) ? newState : arc.weight;
		terminalState = (arc.weight > terminalState) ? arc.weight : terminalState;
	}
	terminalNode.state2 = terminalState;
}

/**
 *
 * @param cut 
 */
void DD::applyOptimalityCutRelaxed(const Cut &cut) {

	nodes[0].state2 = cut.RHS;

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		const auto netArcId = networkPtr->processingOrder[layer-1].second;
		const auto i_NetNodeId = networkPtr->networkArcs[netArcId].tailId;
		const auto q_NetNodeId = networkPtr->networkArcs[netArcId].headId;

		vulint nodesToRemoved;
		vulint nodesToAdded;

		for (const auto id: tree[layer]) {

			auto& node = nodes[id];
			node.state2 = 0; // INFO setting state to zero before applying cut on the node.

			if (node.incomingArcs.size() > 1) {
				for (size_t i = 1; i < nodes[id].incomingArcs.size(); i++) {
					auto inArcId = node.incomingArcs[i];
					auto& arc = arcs[inArcId];
					const auto& parentNode = nodes[arc.tail];

					auto j_NetNodeId = (arc.decision != -1) ? networkPtr->networkArcs[arc.decision].headId : 0;
					arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
					// create new node for this arc.
					DDNode newNode{number.getNext()};
					newNode.state2 = arc.weight + parentNode.state2;
					arc.head = newNode.id;
					newNode.incomingArcs = {arc.id};
					if (networkPtr->hasStateChanged[layer]) newNode.states = networkPtr->stateUpdateMap.at(layer);
					else {
						newNode.states = parentNode.states;
						if (arc.decision != -1) newNode.states.erase(arc.decision);
					}
					// LATER: loop upto n-2 layers
					// process outgoing arcs.
					if (layer != tree.size()-2) {
						newNode.outgoingArcs = {};
						for (auto outArcId : node.outgoingArcs){
							auto outArc = arcs[outArcId];
							if (newNode.states.count(outArc.decision)){
								auto& childNode = nodes[outArc.head];
								DDArc newArc {number.getNext(), newNode.id, childNode.id, outArc.decision};
								newNode.outgoingArcs.push_back(newArc.id);
								childNode.incomingArcs.push_back(newArc.id);
								arcs.insert(make_pair(newArc.id, newArc));
							}
						}
					}
					else {
						newNode.outgoingArcs = {};
						auto outArcId = node.outgoingArcs[0];
						auto& outArc = arcs[outArcId];
						auto& childNode = nodes[outArc.head];
						DDArc newArc {number.getNext(), newNode.id, childNode.id, outArc.decision};
						newArc.weight = outArc.weight;
						arcs.insert(make_pair(newArc.id, newArc));
						childNode.incomingArcs.push_back(newArc.id);
						newNode.outgoingArcs.push_back(newArc.id);
					}
					// insert to list of nodes.
					nodesToAdded.push_back(newNode.id);
					nodes.insert(make_pair(newNode.id, newNode));
				}
			}
			// process first arc.
			auto inArcId = node.incomingArcs[0];
			auto& arc = arcs[inArcId];
			const auto& parentNode = nodes[arc.tail];
			node.incomingArcs = {inArcId}; // reset the incoming arcs.

			if (networkPtr->hasStateChanged[layer]) node.states = networkPtr->stateUpdateMap.at(layer);
			else {
				node.states = parentNode.states;
				if (arc.decision != -1) node.states.erase(arc.decision);
			}
			auto j_NetNodeId = (arc.decision != -1 ) ? networkPtr->networkArcs[arc.decision].headId  : 0;
			arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
			node.state2 = arc.weight + parentNode.state2;
			// remove outArcs that are infeasible.
			if (layer != tree.size()-2){
				auto outgoingArcs = node.outgoingArcs;
				for (auto outArcId : outgoingArcs){
					auto& outArc = arcs[outArcId];
					if (!node.states.count(arcs[outArcId].decision)) {
						// remove that arc
						deleteArcById(outArcId);
					}
				}
			}
		}

		for (const auto id: nodesToRemoved) removeNode(id); // remove nodes marked for deletion
		tree[layer].insert(tree[layer].end(), nodesToAdded.begin(), nodesToAdded.end()); // add node ids to current layer.
	}

	auto& terminalNode = nodes[tree[tree.size()-1][0]];
	double terminalState = std::numeric_limits<double>::lowest();

	for (const auto arcId : terminalNode.incomingArcs){
		auto& arc = arcs[arcId];
		auto newState = nodes[arc.tail].state2;
		arc.weight = (newState < arc.weight) ? newState : arc.weight;
		terminalState = (arc.weight > terminalState) ? arc.weight : terminalState;
	}
	terminalNode.state2 = terminalState;
}

/**
 * Applies optimality cut to the restricted DD.
 * @param cut
 */
double DD::applyOptimalityCutRestrictedLatest(const Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for the sub root node.
	auto& rootNode = nodes[0];
	auto justifiedRHS = cut.RHS;
	for (size_t i = 0; i < rootNode.solutionVector.size(); i++) {
		auto decision = rootNode.solutionVector[i];
		if (decision != -1) {
			auto netArcId = processingOrder[i].second;
			uint iNetId = netArcs[netArcId].tailId;
			uint qNetId = netArcs[netArcId].headId;
			uint jNetId = netArcs[decision].headId;
			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
		}
	}
	rootNode.state2 = justifiedRHS;
#ifdef DEBUG
	cout << "Optimality Actual" << endl;
	// cout << "Optimality cut. Justified RHS: " << justifiedRHS;
#endif

	for (size_t layer = 1; layer < tree.size()-1; layer++) {

		auto netArcId = processingOrder[startTree+layer-1].second;
		auto iNetNodeId = netArcs[netArcId].tailId;
		auto qNetNodeId = netArcs[netArcId].headId;

		for (auto nodeId: tree[layer]) {
			assert(nodes.count(nodeId)>0);
			auto& node = nodes[nodeId];
			auto inArcId = node.incomingArcs[0];
			assert(arcs.count(inArcId)>0);
			auto& inArc = arcs[inArcId];
			const auto& parentNode = nodes[inArc.tail];

			if (parentNode.states.count(inArc.decision)) { //  this condition is unnecesary. since tree is restricted.
				// jNetNodeId
				auto jNetNodeId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
				inArc.weight = (inArc.decision != -1) ? cut.cutCoeff.at(make_tuple(iNetNodeId, qNetNodeId, jNetNodeId))  :0;
				node.state2 = inArc.weight + parentNode.state2;
			}
		}
	}
	// terminal layer
	assert(nodes.count(tree[tree.size()-1][0])> 0);
	auto& terminalNode = nodes[tree[tree.size()-1][0]];
	terminalNode.state2 = 0;
	double terminalState = INT32_MIN;
	for (auto inArcId : terminalNode.incomingArcs) {
		auto& arc = arcs[inArcId];
		const auto& parentNode = nodes[arc.tail];

		// if (arc.weight > parentNode.state2) {
		// 	arc.weight = parentNode.state2;
		// }
		// if(terminalNode.state2 < arc.weight)
		// 	terminalNode.state2 = arc.weight;
		arc.weight = min(arc.weight, parentNode.state2);
		terminalState = max(terminalState, arc.weight);

	}
	// for (auto inArcId : terminalNode.incomingArcs) { // fuse this loop with above.
	// 	if (terminalNode.state2 < arcs[inArcId].weight) {
	// 		terminalNode.state2 = arcs[inArcId].weight;
	// 	}
	// }

#ifdef DEBUG
	cout << "\t terminal state: " << terminalNode.state2 << endl;
#endif
	terminalNode.state2 = terminalState;
	return terminalNode.state2;
}

/**
 * Applies feasibility cut on the restricted DD. Removes infeasible solutions from DD.
 * @param cut
 */
bool DD::applyFeasibilityCutRestrictedLatest(const Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for the sub root node.
	auto& rootNode = nodes[0];
	auto justifiedRHS = cut.RHS;
	for (size_t i = 0; i < rootNode.solutionVector.size(); i++) {
		auto decision = rootNode.solutionVector[i];
		if (decision != -1) {
			auto netArcId = processingOrder[i].second;
			uint iNetId = netArcs[netArcId].tailId;
			uint qNetId = netArcs[netArcId].headId;
			uint jNetId = netArcs[decision].headId;
			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
		}
	}
	rootNode.state2 = justifiedRHS;
#ifdef DEBUG
	 cout << "Feasibility cut. Acutal" << endl;
#endif

	// auto lowerBounds = helperFunction(network, cut);
	// lowerBounds.push_back(0);

	vulint nodesToRemove;
	// nodesToRemove.reserve(MAX_WIDTH);

	for (size_t layer = 1; layer < tree.size()-1; layer++) {

		auto netArcId = processingOrder[startTree+layer-1].second;
		auto iNetNodeId = netArcs[netArcId].tailId;
		auto qNetNodeId = netArcs[netArcId].headId;

		nodesToRemove.clear(); // not using lower bound heuristic
		auto globalPosition = startTree + layer -1;

		for(auto nodeId : tree[layer]) {
			assert(nodes.count(nodeId)>0);
			auto& node = nodes[nodeId];
			auto inArcId = node.incomingArcs[0];
			assert(arcs.count(inArcId)>0);
			auto& inArc = arcs[inArcId];
			const auto& parentNode = nodes[inArc.tail];

			if (inArc.decision != -1) {
				auto jNetNodeId = netArcs[inArc.decision].headId;
				inArc.weight = cut.cutCoeff.at(make_tuple(iNetNodeId, qNetNodeId, jNetNodeId));
				node.state2 = parentNode.state2 + inArc.weight;
				// if ((newState + lowerBounds[node.globalLayer]) >= 0)
				// 	node.state2 = newState;
				// else {
				// 	node.state2 = newState;
				// 	nodesToRemove.push_back(nodeId);
				// }
			}
			else {
				inArc.weight = 0;
				node.state2 = parentNode.state2;
				// if ((node.state2 + lowerBounds[node.globalLayer]) < 0)
				// 	nodesToRemove.push_back(nodeId);
			}
			if(node.state2 < -0.5) nodesToRemove.push_back(node.id);
		}

		// if (!nodesToRemove.empty()) {
		// 	#ifdef DEBUG
		// 		cout << " # of nodes that are infeasible in layer "<< layer << " : " << nodesToRemove.size() << endl;
		// 	#endif
		// 	if (nodesToRemove.size() == tree[layer].size()) {
		// 		#ifdef DEBUG
		// 			cout << "Removing entire layer.. " << layer << endl;
		// 		#endif
		// 		// set isFeasible flag to false
		// 		isInFeasible = true;
		// 		// TODO do not remove entire tree/nodes, set flags and discard current semi root.
		// 	}
		// 	batchRemoveNodes(nodesToRemove);
		// }
	}

	if (nodesToRemove.size() == tree[tree.size()-2].size()) {
		// removing entire layer removes entire tree. update flags instead
		#ifdef DEBUG
		cout << "Entire layer is deleted. DD is infeasible." << endl;
		#endif
		isInFeasible = true;
		return false;
	}
#ifdef DEBUG
	if (!nodesToRemove.empty()) cout << " Removing " << nodesToRemove.size() << " nodes from the tree" << endl;
#endif
	// cout << endl;
	if (!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove);

	return true;
}


double DD::applyOptimalityCutHeuristic(const Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for the sub root node.
	auto& rootNode = nodes[0];
	auto justifiedRHS = cut.RHS;
	for (size_t i = 0; i < rootNode.solutionVector.size(); i++) {
		auto decision = rootNode.solutionVector[i];
		if (decision != -1) {
			auto netArcId = processingOrder[i].second;
			uint iNetId = netArcs[netArcId].tailId;
			uint qNetId = netArcs[netArcId].headId;
			uint jNetId = netArcs[decision].headId;
			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
		}
	}
	rootNode.state2 = justifiedRHS;
#ifdef DEBUG
	cout << "Optimality cut heuristic" << endl;
#endif

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArcId = processingOrder[startTree+layer-1].second;
		auto iNetNodeId = netArcs[netArcId].tailId;
		auto qNetNodeId = netArcs[netArcId].headId;

		for (auto nodeId : tree[layer]) {
			assert(nodes.count(nodeId)>0);
			auto& node = nodes[nodeId];
			double newState = std::numeric_limits<double>::lowest();

			for (auto inArcId : node.incomingArcs) {
				assert(arcs.count(inArcId)> 0);
				auto& arc = arcs[inArcId];
				assert(nodes.count(arc.tail)>0);
				auto& parentNode = nodes[arc.tail];
				auto jNetNodeId = (arc.decision != -1) ? netArcs[arc.decision].headId : 0;
				arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(iNetNodeId, qNetNodeId, jNetNodeId))  :0;
				if (newState <= (arc.weight + parentNode.state2)) {
					newState = arc.weight + parentNode.state2;
				}
			}
			node.state2 = newState;
		}
	}

	double terminalState = std::numeric_limits<double>::lowest();

	auto& terminalNode = nodes[tree[tree.size()-1][0]];
	for (auto inArcId : terminalNode.incomingArcs) {
		auto& inArc = arcs[inArcId];
		auto& parentNode = nodes[inArc.tail];
		// if (inArc.weight > parentNode.state2) {
		// 	inArc.weight = parentNode.state2;
		// }
		// min of terminal arc weight and parent state. and maximum of terminal arcs weights.
		inArc.weight = min(inArc.weight, parentNode.state2);
		terminalState = std::max(terminalState, inArc.weight);
	}
	// for (auto inArcId : terminalNode.incomingArcs) {
	// 	if (terminalState < arcs[inArcId].weight) {
	// 		terminalState = arcs[inArcId].weight;
	// 	}
	// }
	terminalNode.state2 = terminalState;
#ifdef DEBUG
	cout << "\t terminal state: " << terminalState << endl;
#endif
	// cout << "Applied Optimality cut heuristically" << endl;
	return terminalState;
}

/**
 * Applies given feasiblity cut on  the DD. Returns boolean true if the
 * DD remains feasible after applying the cut, boolean false otherwise.
 * This function will not prune solutions and should only be applied on the relaxed tree.
 *
 * @param cut
 * @return
 */
bool DD::applyFeasibilityCutHeuristic(const Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for the sub root node.
	auto& rootNode = nodes[0];
	auto justifiedRHS = cut.RHS;
	for (size_t i = 0; i < rootNode.solutionVector.size(); i++) {
		auto decision = rootNode.solutionVector[i];
		if (decision != -1) {
			auto netArcId = processingOrder[i].second;
			uint iNetId = netArcs[netArcId].tailId;
			uint qNetId = netArcs[netArcId].headId;
			uint jNetId = netArcs[decision].headId;
			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
		}
	}

	rootNode.state2 = justifiedRHS;

	for (size_t layer = 1; layer < tree.size() -1; layer++) {
		auto theArc = processingOrder[startTree+layer-1].second;
		const auto& iNetNode = netArcs[theArc].tailId;
		const auto& qNetNode = netArcs[theArc].headId;

		for (auto nodeId : tree[layer]) {
			double newState = std::numeric_limits<double>::lowest();
			auto& node = nodes[nodeId];
			for (auto inArcId : node.incomingArcs) {
				auto& inArc = arcs[inArcId];
				const auto& parentNode = nodes[inArc.tail];
				const auto& jNetNode = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
				inArc.weight = (inArc.decision != -1) ? cut.get(iNetNode, qNetNode, jNetNode) : 0;
				newState = (inArc.weight + parentNode.state2 > newState) ? inArc.weight + parentNode.state2 : newState;
			}
			node.state2 = newState;
		}
	}

	vector<ulint> nodesToRemove; // remove nodes that have negative state.
	for (auto nodeId : tree[tree.size()-2]) {
		if (nodes[nodeId].state2 < 0) nodesToRemove.push_back(nodeId);
	}
	if (nodesToRemove.size() == tree[tree.size()-2].size()) return false;
	if (!nodesToRemove.empty())batchRemoveNodes(nodesToRemove);
	return true;
	// terminal layer
	vui arcToDeleted;
	double state = std::numeric_limits<double>::lowest();
	auto terminalNode = nodes[tree[tree.size()-1][0]];
	for (auto inArcId : terminalNode.incomingArcs) {
		auto& arc = arcs[inArcId];
		const auto& parentNode = nodes[arc.tail];
		// arc.weight = (arc.weight > parentNode.state2) ? parentNode.state2 : arc.weight;
		// state = (state < arc.weight) ? arc.weight : state;
		arc.weight = min(arc.weight,parentNode.state2);
		state = max(state, arc.weight);
		if (arc.weight < 0) arcToDeleted.push_back(inArcId);
	}

	if (state < 0) return false;
	cout << "Feasibility Heuristic: deleting..." << arcToDeleted.size() << " arcs from terminal layer" << endl;
	for (auto arcId : arcToDeleted) { deleteArcById(arcId);} // should work.

#ifdef DEBUG
	cout << "\t terminal state: " << state << endl;
#endif
	terminalNode.state2 = state;
	return (state >= 0);
}


/**
 * Compiles the decision diagram tree with given root as the root.
 * @param root - (sub-) root of the decision diagram.
 * @return - vector of nodes that are in exact cutset if tree becomes Restricted.
 *
 * NOTE: this function will return the cutset once. If the tree is exact,
 * function returns empty, thus the return parameter "optional".
 */
// optional<vector<Node_t>> RestrictedDD::compile(DDNode root) {
// 	// compile tree based on the root.
//
// 	const auto& processingOrder = networkPtr->processingOrder;
// 	const auto& stateUpdateMap = networkPtr->stateUpdateMap;
//
// 	vui currentLayer;
//
// 	startTree = root.globalLayer;
// 	// build root node and insert it to the layer.
// 	root.incomingArcs.clear();
// 	root.outgoingArcs.clear();
// 	root.id = 0;
// 	root.nodeLayer = 0;
// 	uint exactLayer = 0;
//
// 	nodes.insert(make_pair(0, root));
// 	currentLayer.push_back(0);
//
// 	tree.push_back(currentLayer);
//
// 	auto start = processingOrder.begin() + startTree;
// 	auto end = processingOrder.end();
//
// 	uint index = 0;
// 	uint nextSize = 1;
//
// 	for (; start != end; ++start) {
// 		auto[a,b] = *start;
//
// 		if (stateUpdateMap.count(a)) {
// 			const auto& newStates = stateUpdateMap.at(a);
// 			updateStates(currentLayer, newStates);
// 		}
// 		bool isExact = true;
// 		vui nextLayer = buildNextLayer(currentLayer, networkPtr->hasStateChanged[a+1], isExact, nextSize);
// 		if (isExact) exactLayer++;
// 		++index;
// 		tree.push_back(nextLayer);
// 		currentLayer = std::move(nextLayer);
// 	}
//
// 	// terminal layer.
// 	vui terminalLayer;
// 	DDNode terminalNode{number++};
// 	terminalNode.nodeLayer = ++index;
//
// 	// current layer points to pre-terminal layer.
// 	for (const auto& id : currentLayer) {
// 		uint arcId = number++;
// 		auto& parentNode = nodes[id];
// 		DDArc arc{arcId, id, terminalNode.id, 1};
// 		arc.weight = std::numeric_limits<double>::max();
// 		terminalNode.incomingArcs.push_back(arcId);
// 		parentNode.outgoingArcs.push_back(arcId);
// 		arcs.insert(make_pair(arcId, arc));
// 	}
//
// 	nodes.insert(make_pair(terminalNode.id, terminalNode));
// 	terminalLayer.push_back(terminalNode.id);
// 	tree.push_back(terminalLayer);
//
// #ifdef DEBUG
// 	displayStats();
// #endif
// 	if (exactLayer == index-1) return {};
// 	return generateExactCutSet(exactLayer);
// }
//
// /**
//  * Builds next layer given the current layer.
//  * @param currentLayer
//  * @param hasStateChanged
//  * @param isExact
//  * @param nextSize
//  * @return
//  */
// vui RestrictedDD::buildNextLayer(const vui &currentLayer, bool hasStateChanged, bool &isExact, uint& nextSize) {
// 	vui nextLayer;
// 	// reserve space based on the nextSize.
//
// 	return nextLayer;
//
// }
//
//
// /**
//  * Computes and returns a vector of nodes that are in the cutset.
//  * @param exactLayer - Id of exact layer in the tree.
//  * @return vector of nodes that are in exact cutset.
//  */
// vector<Node_t> RestrictedDD::generateExactCutSet(uint exactLayer) const {
// 	vector<Node_t> cutsetNodes;
// 	cutsetNodes.reserve(tree[exactLayer].size());
// 	for (const auto id: tree[exactLayer]) {
// 		const auto& node = nodes.at(id);
// 		vi states{node.states.begin(), node.states.end()};
// 		cutsetNodes.emplace_back(states, computePathForNode(id),
// 			std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max(), node.globalLayer);
// 	}
// 	return cutsetNodes;
// }
//
// /**
//  * Returns solution vector for the given node in the eact layer. this function should only
//  * be used for a node in the exact cutset.
//  * @param nodeId - id of the node in the cutset.
//  * @return - vector of decisions from root to the given node.
//  */
// vi RestrictedDD::computePathForNode(uint nodeId) const {
// 	DDNode node = nodes.at(nodeId);
// 	vi solution;
//
// 	while (!node.incomingArcs.empty()) {
// 		// TODO: optimize without creating new node each time.
// 	}
//
// 	return solution;
// }
//
// /**
//  * Removes the node and its associated nodes from the nodes map and updates the
//  * tree vector if deletion is not occurring in batch.
//  * @param id Id of the node to be removed.
//  * @param isBatch updates the tree vector if parameter is false. Default value is false.
//  */
// void RestrictedDD::removeNode(uint id, bool isBatch){
// 	// /*
// 	//  * NOTE:
// 	//  */
// 	// // if current node has multiple parents and multiple children, do this.
// 	// auto& node = nodes[id];
// 	// auto outArcs = node.outgoingArcs;
// 	// assert(!node.outgoingArcs.empty());
// 	// for (auto childArcId : outArcs){
// 	// 	assert(arcs.count(childArcId));
// 	// 	auto& childArc = arcs[childArcId];
// 	// 	assert(nodes.count(childArc.head));
// 	// 	auto& child = nodes[childArc.head];
// 	// 	deleteArc(node, childArc, child);
// 	// 	if (child.incomingArcs.empty()) topDownDelete(child.id); // orphan node
// 	// }
// 	// //if (node.incomingArcs.size() > 1) { // multiple parents
// 	// auto incomingArcs = node.incomingArcs;
// 	// assert(!node.incomingArcs.empty());
// 	// for (auto arcId: incomingArcs){
// 	// 	assert(arcs.count(arcId));
// 	// 	auto& arc = arcs[arcId];
// 	// 	assert(nodes.count(arc.tail));
// 	// 	auto& parentNode = nodes[arc.tail];
// 	// 	deleteArc(parentNode, arc, node);
// 	// 	if (parentNode.outgoingArcs.empty()) bottomUpDelete(parentNode.id); // orphan node
// 	// }
// 	// assert(node.incomingArcs.empty() && node.outgoingArcs.empty());
// 	// // actual delete.
// 	// deleteNode(node);
// 	// assert(!nodes.count(id));
// 	//
// 	// if (!isBatch) updateTree();	// no batch deletion.
// }
//
//
// void RestrictedDD::batchRemoveNodes(const vui& nodeIds) {
//
// 	// for (auto id: nodeIds) removeNode(id, true);
//
// 	// updateTree();
// }
//
// double RestrictedDD::applyOptimalityCut(const Cut &cut) noexcept {
//
// 	const auto& processingOrder = networkPtr->processingOrder;
// 	const auto& netArcs = networkPtr->networkArcs;
// 	// compute justified RHS for (sub)-root.
// 	auto& rootNode = nodes[0];
// 	auto justifiedRHS = cut.RHS;
// 	for (size_t i = 0; i < rootNode.solutionVector.size(); i++) {
// 		auto decision = rootNode.solutionVector[i];
// 		if (decision != -1) {
// 			auto netArcId = processingOrder[i].second;
// 			uint iNetId = netArcs[netArcId].tailId;
// 			uint qNetId = netArcs[netArcId].headId;
// 			uint jNetId = netArcs[decision].headId;
// 			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
// 		}
// 	}
// 	rootNode.state2 = justifiedRHS;
//
// 	for (size_t layer = 1; layer < tree.size()-1; layer++) {
// 		const auto netArcId = processingOrder[startTree+layer-1].second;
// 		const auto iNetId = netArcs[netArcId].tailId;
// 		const auto qNetId = netArcs[netArcId].headId;
// 		for (auto nodeId : tree[layer]) {
// 			auto& node = nodes[nodeId];
// 			auto inArcId = node.incomingArcs[0];
// 			auto& inArc = arcs[inArcId];
// 			const auto& parentNode = nodes[inArc.tail];
// 			// removed a condition to check this arc's decision in parent node's state.
// 			// since this is restricted DD, decision must be in parent's state.
// 			assert(parentNode.states.count(inArc.decision)>0);
// 			auto jNetId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
// 			inArc.weight = (inArc.decision != -1) ? cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId)) : 0;
// 			node.state2 = inArc.weight + parentNode.state2;
// 		}
// 	}
// 	// terminal layer.
// 	auto& terminalNode = nodes[tree.back().back()];
// 	terminalNode.state2 = 0;
// 	for (auto inArcId : terminalNode.incomingArcs) {
// 		auto& arc = arcs[inArcId];
// 		const auto& parentNode = nodes[arc.tail];
// 		// if (arc.weight > parentNode.state2) arc.weight = parentNode.state2;
// 		arc.weight = (arc.weight > parentNode.state2) ? parentNode.state2 : arc.weight;
// 	}
// 	// LATER: fuse below loop with above.
// 	for (auto inArcId : terminalNode.incomingArcs) {
// 		terminalNode.state2 = (terminalNode.state2 < arcs[inArcId].weight) ? arcs[inArcId].weight : terminalNode.state2;
// 	}
// 	return terminalNode.state2;
// }
//
// bool RestrictedDD::applyFeasibilityCut(const Cut& cut) noexcept {
//
// 	const auto& processingOrder = networkPtr->processingOrder;
// 	const auto& netArcs = networkPtr->networkArcs;
//
// 	// compute justified RHS for (sub-) root.
// 	auto& rootNode = nodes[0];
// 	auto justifiedRHS = cut.RHS;
// 	for(size_t i = 0; i < rootNode.solutionVector.size(); i++) {
// 		auto decision = rootNode.solutionVector[i];
// 		if (decision != -1) {
// 			auto netArcId = processingOrder[i].second;
// 			uint iNetId = netArcs[netArcId].tailId;
// 			uint qNetId = netArcs[netArcId].headId;
// 			uint jNetId = netArcs[decision].headId;
// 			justifiedRHS += cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
// 		}
// 	}
// 	rootNode.state2 = justifiedRHS;
//
// 	vui nodesToRemove;
//
// 	for (size_t layer = 1; layer < tree.size()-1; layer++) {
// 		auto netArcId = processingOrder[startTree+layer-1].second;
// 		auto iNetId = netArcs[netArcId].tailId;
// 		auto qNetId = netArcs[netArcId].headId;
//
// 		nodesToRemove.clear();
// 		auto globalPosition = startTree + layer-1;
// 		for (auto nodeId : tree[layer]) {
// 			auto& node = nodes[nodeId];
// 			auto inArcId = node.incomingArcs[0];
// 			auto& inArc = arcs[inArcId];
// 			const auto& parentNode = nodes[inArc.tail];
//
// 			if (inArc.decision != -1) {
// 				auto jNetId = netArcs[inArc.decision].headId;
// 				inArc.weight = cut.cutCoeff.at(make_tuple(iNetId, qNetId, jNetId));
// 				node.state2 = parentNode.state2 + inArc.weight;
// 			}
// 			else {
// 				inArc.weight = 0;
// 				node.state2 = parentNode.state2;
// 			}
// 			if (node.state2 < -0.5) nodesToRemove.push_back(nodeId);
// 		}
// 	}
//
// 	if (nodesToRemove.size() == tree[tree.size()-2].size()) {
// 		// entire layer is to be removed -> tree is infeasible.
// #ifdef DEBUG
// 		cout << "Entire last layer is deleted. DD is infeasible." << endl;
// #endif
// 		isTreeDeleted = true;
// 		return false;
// 	}
//
// 	if(!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove);
//
// 	return true;
// }

/************************************** BEGIN: RESTRICTED TREE ************************************/

optional<vector<Inavap::Node>> Inavap::RestrictedDD::buildTree(Inavap::Node root) {

	// node parameter should be initialized before calling this function. node should contain states.
	// the given node parameter will be inserted as the root of the tree.

	const auto& arcOrder = networkPtr->processingOrder;
	const auto& stateUpdateMap = networkPtr->stateUpdateMap;

	// WIDTH = WIDTH_;

	// set the node to the root.
	vector<uint> currentLayer; // should be root layer.
	startTree = root.globalLayer;
	rootSolution = root.solutionVector;

	RDDNode node{root};
	currentLayer.push_back(node.id);
	// insert root layer to tree.
	tree.push_back(currentLayer);
	nodes.insert(std::make_pair(node.id, node));

	auto start = arcOrder.begin() + startTree;
	auto end = arcOrder.end();

	tree.reserve(std::distance(start, end)); // pre-allocate in vector.

	uint index = 0; // for node layer.
	uint8_t isExact = true;
	uint nextLayerSize = 0; // stores size of the next layer.
	uint exactLayer = 0; // index of exact layer.

	// build tree until terminal layer.
	for (; start != end; ++start){

		auto[a,b] = *start;
		if (stateUpdateMap.contains(a)) {
			// update state of each node in the layer if state changes next.
			const auto& newStates = stateUpdateMap.at(a);
			vector<int16_t> statesVector;
			for (auto state: newStates) {
				statesVector.push_back(static_cast<int16_t>(state));
			}
			// const vector<shi> statesVector{newStates.begin(), newStates.end()};
			updateStates(currentLayer, statesVector);
		}

		auto stateChangesNext = networkPtr->hasStateChanged[a+1];
		vector<uint> nextLayer = buildNextLayer(currentLayer,nextLayerSize,
			stateChangesNext, isExact);
		if (isExact) exactLayer++; // at last, this number should be exact layer number.
		++index;
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}

	terminalId = ++lastInserted;
	RDDNode terminalNode {terminalId};
	terminalNode.nodeLayer = ++index;
	// terminal node layer.
	vector<uint> terminalLayer;
	terminalLayer.push_back(terminalId);
	tree.push_back(terminalLayer);

	/* NOTE:: current layer now points to last layer of tree, and nextLayerSize is the number of arcs created
	 * in terminal layer. Terminal node doesn't store the incoming arcs, instead terminal arcs are stored and
	 * modified separately in the RestrictedDD class.
	 */
	terminalInArcs.reserve(nextLayerSize);
	for (auto id : currentLayer) {
		// create new arc for each node in current layer. we do not create separate arc for '0'.
		auto& parentNode = nodes[id];
		uint arcId = ++lastInserted; // arc and the head node have different Id.
		DDArc arc{arcId, id, terminalId, 1};
		arc.weight = std::numeric_limits<double>::max();
		parentNode.outgoingArcs.push_back(arcId);
		terminalInArcs.push_back(arcId);
		arcs.insert(make_pair(arcId, arc));
	}

	nodes.insert(make_pair(terminalId, terminalNode));

	// if (exactLayer == index-1) return {};
	if (isExact) return {};
	status = 1;
	return generateExactCutSet(exactLayer);
}

/**
 * Returns the path from global root to the given node.
 */
vector<int16_t> Inavap::RestrictedDD::getSolutionForNode(uint id) const {

	// since every node has one parent, simply traverse from current node till root.

	vector<int16_t> reversePath;
	const RDDNode *current = &nodes.at(id);
	while (current->nodeLayer) { // iterate recursively until root.
		const auto& inArc = arcs.at(current->incomingArc);
		reversePath.push_back(inArc.decision);
		current = &nodes.at(inArc.tail);
	}
	// add root's solution.
	vector result(rootSolution);
	result.insert(result.end(), reversePath.rbegin(), reversePath.rend());
	return result;
}

/**
 * Returns the vector of cutset nodes from the restricted tree.
 * @param layer
 * @return
 */
vector<Inavap::Node> Inavap::RestrictedDD::generateExactCutSet(uint layer) const{
	// find exact layer first, then generate solution for each node in cutset.
	vector<Node> cutsetNodes;
	cutsetNodes.reserve(tree[layer].size());

	for (const auto id: tree[layer]) {
		const auto& node = nodes.at(id);
		cutsetNodes.emplace_back(node.states, getSolutionForNode(id),
			std::numeric_limits<double>::lowest(),std::numeric_limits<double>::max(), node.globalLayer
		);
	}
	return cutsetNodes;
}

vector<int16_t> Inavap::RestrictedDD::getSolution() const noexcept {
	// find the max node and its id from terminal arcs.
	double maxVal = std::numeric_limits<double>::lowest();
	uint maxNodeId = 0;

	for (const auto arcId : terminalInArcs) {
		const auto& arc = arcs.at(arcId);
		if (arc.weight > maxVal) {
			maxVal = arc.weight;
			maxNodeId = arc.tail;
		}
	}
	return getSolutionForNode(maxNodeId);
	// vector<int16_t> reversePath;
	// // reversePath.reserve();
	//
	// const RDDNode *current = &nodes.at(maxNodeId);
	//
	// while (current->nodeLayer) {
	// 	const auto& inArc = arcs.at(current->incomingArc);
	// 	reversePath.push_back(inArc.decision);
	// 	current = &nodes.at(inArc.tail);
	// }
}

vector<uint> Inavap::RestrictedDD::buildRestrictedLayer(const vector<uint> &currentLayer) {
	// tree is not exact.
	// strategy: for each node create one child. skip -1 for a node, unless otherwise.
	vector<uint> nextLayer;
	nextLayer.reserve(WIDTH);

	for (const auto id: currentLayer) {
		// invariant: node Id and the arc incoming to the node will have same id.
		auto& parent = nodes[id];
		int16_t decision = parent.states.back(); // get last decision from parent states.
		auto childStates = parent.states;
		if (decision != -1) { // remove decision from child states.
			childStates.erase(find(childStates.begin(), childStates.end(), decision));
		}
		uint index = ++lastInserted;
		// create new node with this decision.
		DDArc arc{index, id, index, decision};
		parent.outgoingArcs.push_back(arc.id);
		RDDNode childNode {index};
		childNode.globalLayer = parent.globalLayer + 1;
		childNode.nodeLayer = parent.nodeLayer + 1;
		childNode.incomingArc = index;
		childNode.states = std::move(childStates);
		nextLayer.push_back(childNode.id);
		nodes.insert(make_pair(index, childNode));
		arcs.insert(make_pair(index, arc));
	}
	return nextLayer;
}

vector<uint> Inavap::RestrictedDD::buildNextLayer(const vector<uint> &currentLayer,
                                                  uint& nextLayerSize, uint8_t stateChangesNext, uint8_t &isExact) {

	// build the next layer until MAX WIDTH is reached.

	// isExact = true;
	vector<uint> nextLayer;
	nextLayer.reserve(min(WIDTH, nextLayerSize));
	nextLayerSize = 0;

	if (isExact) {
		uint count = 0;
		// tree is still exact, build complete layer.
		for (auto id : currentLayer) {
			RDDNode &parentNode = nodes[id];
			const auto parentStates = parentNode.states;

			// states contain -1, do not prioritize it unless needed. assuming states are in ascending order (typically).
			for (auto start = parentNode.states.rbegin(); start != parentNode.states.rend(); ++start) {
				auto decision = *start;
				if (count >= WIDTH) { isExact = false; return nextLayer;}
				auto index = ++lastInserted;
				RDDNode node{index};
				DDArc arc{index, parentNode.id, node.id, decision};
				parentNode.outgoingArcs.push_back(arc.id);
				node.states = parentStates;
				if (decision != -1) node.states.erase(find(node.states.begin(), node.states.end(), decision));
				node.incomingArc = arc.id;
				node.nodeLayer = parentNode.nodeLayer+1;
				node.globalLayer = parentNode.globalLayer + 1;
				nextLayerSize += node.states.size(); // update next layer size.
				// insert node and arc to containers.
				nodes.insert(make_pair(index, node));
				arcs.insert(make_pair(index, arc));
				count++;
				nextLayer.push_back(index);
			}
		}
		return nextLayer;
	}

	// tree is not exact, restricted layer.
	// strategy. do not build

	return buildRestrictedLayer(currentLayer);
}

[[always_inline]] void Inavap::RestrictedDD::updateStates(const vector<uint> &currentLayer, const vector<int16_t> &nextLayerStates) {

	for_each(currentLayer.begin(), currentLayer.end(), [this, &nextLayerStates](uint id) {
		auto& node = nodes.at(id);
		node.states = nextLayerStates;
	});
}

void Inavap::RestrictedDD::deleteArc(RDDNode& tailNode, DDArc& arc, RDDNode& headNode) noexcept {
	auto pos = std::find(tailNode.outgoingArcs.begin(), tailNode.outgoingArcs.end(), arc.id);
	if (pos != tailNode.outgoingArcs.end()) tailNode.outgoingArcs.erase(pos);
	// NOTE:  no need to update the incoming of head node. going to be deleted anyway.
	arcs.erase(arc.id);
}

/**
 * Deletes node from the node's container. Calls node container's erase function.
 */
void Inavap::RestrictedDD::deleteNode(RDDNode& node) {
	deletedNodeIds.emplace_back(node.id);
	nodes.erase(node.id);
	numberOfDeletedNodes++;
}

void Inavap::RestrictedDD::updateTree() {

	// sort by index (indirectly sorts by layer).
	sort(deletedNodeIds.begin(), deletedNodeIds.end());
	// uint count = deletedNodeIds.size();

	/* since both the deleted ids and layer is sorted, existence of an deleted Id in the given layer
	 * can be performed in O(1) time : compare with last element of layer.
	 *
	 * Invariant: current max in deleted Ids <= max Id in the layer.
	 */
	auto f = [this](vui& ids) {
		if (this->deletedNodeIds.empty()) return;
		//  vector doesn't support erasing with reverse iterator. create a mask for elements that are removed.
		vector<uint8_t> mask;
		mask.reserve(ids.size());

		for_each(ids.rbegin(), ids.rend(), [this,&mask](uint id) {
			// NOTE: the above invariant must be hold all time in order trick to work.
			if (!deletedNodeIds.empty() && id==this->deletedNodeIds.back()) {
				mask.push_back(1);
				this->deletedNodeIds.pop_back();
				// count--;
			}
			else mask.push_back(0); // edge case handled here.
		});

		std::reverse(mask.begin(), mask.end());
		size_t index = 0;
		ids.erase(std::remove_if(ids.begin(), ids.end(),
			[&mask,&index](uint id) { return static_cast<bool>(mask[index++]); }), ids.end());
	};
	for_each(tree.rbegin(), tree.rend(), f);
	deletedNodeIds.clear();
}

/**
 * Deletes the entire subtree of the given node recursively till the terminal. deletes the given node from the
 * node's container.
 */
void Inavap::RestrictedDD::topDownDelete(RDDNode& node) {
	// not a terminal node.
	auto outArcs = node.outgoingArcs;

	for (auto childArcId: outArcs) {
		auto& childArc = arcs[childArcId];
		auto& child = nodes[childArc.head];
		deleteArc(node, childArc, child);
		if (child.id != terminalId) topDownDelete(child); // do not delete terminal node.
	}
	deleteNode(node);
}

/**
 * Deletes the node recursively until the root. deletes the given node from the node's container.
 */
void Inavap::RestrictedDD::bottomUpDelete(RDDNode& node) {
	// INVARIANT: this node should not contain any outgoing arcs.
	// remove incoming arc and call bottom-up delete recursively.
	auto& arc = arcs[node.incomingArc];
	auto& parent = nodes[arc.tail];
	deleteArc(parent, arc, node);
	deleteNode(node);
	if (parent.outgoingArcs.empty()) bottomUpDelete(parent);
}

/**
 * Deletes the given node and its subtree, and associated parent(s) (if eligible) from the node's container.
 * Also updates the tree in-case of non-batch deletions.
 * @param id - Id of the node to be deleted.
 * @param isBatch - If the function is called in a batch deletions?
 */
void Inavap::RestrictedDD::removeNode(uint id, bool isBatch) {

	// just following Erfan's strategy, just remove the arc between the node and the terminal node.
	auto&  temp_node = nodes[id];
	if (temp_node.outgoingArcs.empty()) return; // already deleted node.
	auto& temp_arc = arcs[temp_node.outgoingArcs[0]];
	auto& terminal_node = nodes[terminalId];
	terminalInArcs.erase(find(terminalInArcs.begin(), terminalInArcs.end(), temp_arc.id));
	deleteArc(temp_node, temp_arc, terminal_node);
	temp_node.outgoingArcs.clear();

	return;
	/* Deletes the incoming arc from parent, update its outgoing arcs and start top-down delete
	 * for the node. Top-down delete should internally delete the given node and update the node
	 * container. Caller is responsible to check the given id is not the root id.
	 */

	auto& node = nodes[id];

	/* NOTE: This function is only called on the nodes from the last layer (before terminal) that have only one incoming
	 * and one outgoing arc. Incoming arc ids of the terminal node is stored separately, thus update terminalInArcs
	 * vector explicitly.
	 */
	auto& outgoingArc = arcs[node.outgoingArcs[0]];
	auto& terminalNode = nodes[terminalId];
	terminalInArcs.erase(find(terminalInArcs.begin(), terminalInArcs.end(), outgoingArc.id));
	deleteArc(node, outgoingArc, terminalNode); // will not update incoming of terminalNode.

	// delete parent arc. must contain only one parent arc except the root.
	// check for root node.
	auto& incomingArc = arcs[node.incomingArc];
	auto& parentNode = nodes[incomingArc.tail];
	deleteArc(parentNode, incomingArc, node);

	/* According to the way that feasibility cut is applied, we only remove nodes from the last layer.
	 * This will not have any performance benefit if we call top-down delete on these nodes.
	 *
	 * Uncomment below line if new strategy in applying feasibility cut and comment the deleteNode() line.
	 */
	// topDownDelete(node); // delete subtree completely from the container.
	deleteNode(node);
	// delete parent(s) if no children.
	if (parentNode.outgoingArcs.empty()) bottomUpDelete(parentNode);

	if (!isBatch) updateTree(); // batch delete?
}

void Inavap::RestrictedDD::batchRemoveNodes(vector<uint> &ids) {
	//
	for (auto id: ids) removeNode(id, true); // remove nodes without updating tree.
	updateTree();
}

// TODO: store Q's offset from previous call and use optimized search strategy by passing Q's offset in the key.

/**
 * Applies given optimality cut to the tree and returns the lower bound for the tree.
 * @param cut - Optimality cut to apply.
 * @return - lower bound for the tree.
 */
double Inavap::RestrictedDD::applyOptimalityCut(const Inavap::Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& rootNode = nodes[0];
	size_t i = 0;
	rootNode.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](double val , int16_t decision) { // pass reference to val?
			if (decision == -1) {i++;return val;}
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			auto key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);

	// apply function to each layer.

	size_t index = i; // start from index where the root solution ends.

	auto processLayer = [this, &processingOrder, &netArcs, &cut, &index](const auto& layer) {
		const auto netArcId = processingOrder[index++].second;
		const auto iNetId = netArcs[netArcId].tailId;
		const auto qNetId = netArcs[netArcId].headId;

		// process each node in the given layer. update the respective containers.
		auto processNode =  [&cut, iNetId, qNetId, &netArcs, this](auto nodeId) {
			auto& node = nodes.at(nodeId);
			auto& inArc = arcs.at(node.incomingArc);
			const auto& parentNode = nodes.at(inArc.tail);
			// early stopping if decision == -1 without changing the weight to zero.
			if (inArc.decision == -1) {inArc.weight = 0; node.state2 = parentNode.state2; return; }
			auto jNetId = netArcs[inArc.decision].headId;
			auto key = getKey(qNetId, iNetId, jNetId);
			inArc.weight = cut.get(key);
			// if (inArc.decision != -1) {
			// 	auto jNetId = netArcs[inArc.decision].headId;
			// 	auto key = getKey(qNetId, iNetId, jNetId);
			// 	inArc.weight = cut.get(key);
			// }
			// auto jNetId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
			// inArc.weight = (inArc.decision != -1) ? cut.get(iNetId, qNetId, jNetId) : 0;
			node.state2 = inArc.weight + parentNode.state2;
		};
		for_each(layer.begin(), layer.end(), processNode);
	};

	size_t offset = tree.size()-1;
	for_each(tree.begin()+1, tree.begin()+offset, processLayer);

	// terminal layer.
	auto& terminalNode = nodes[terminalId];
	double terminalState = INT32_MIN;

	for_each(terminalInArcs.begin(), terminalInArcs.end(),[this, &terminalState](auto arcId) {
		auto& arc = arcs.at(arcId);
		const auto& parentNode = nodes.at(arc.tail);
		arc.weight = std::min(arc.weight, parentNode.state2);
		terminalState = std::max(terminalState, arc.weight);
	});
	terminalNode.state2 = terminalState;
	return terminalState;
}


uint8_t Inavap::RestrictedDD::applyFeasibilityCut(const Inavap::Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& rootNode = nodes[0];
	size_t i = 0;
	// compute justified RHS for the root node.
	rootNode.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](double val , const int decision) {
			if (decision == -1) {i++; return val;}
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			uint64_t key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);

	// nodes to remove.
	vui nodesToRemove;
	size_t index = i;// startTree; // index of element in processing order vector. TODO verify this.

	auto processLayer = [this, &processingOrder, &netArcs, &cut, &index, &nodesToRemove](const auto& layer) {
		const auto netArcId = processingOrder[index++].second;
		const auto iNetId = netArcs[netArcId].tailId;
		const auto qNetId = netArcs[netArcId].headId;

		nodesToRemove.clear();

		auto processNode =  [&cut, iNetId, qNetId, &netArcs, &nodesToRemove, this](auto nodeId) {
			auto& node = nodes.at(nodeId);
			auto& inArc = arcs.at(node.incomingArc);
			const auto& parentNode = nodes.at(inArc.tail);
			// early stopping if decision == -1 without changing the arc weight.
			if (inArc.decision == -1) {
				inArc.weight = 0;
				node.state2 = parentNode.state2;
				return;
			}
			// inArc.weight = 0;
			// if (inArc.decision != -1) {
			// 	auto jNetId = netArcs[inArc.decision].headId;
			// 	auto key = getKey(qNetId, iNetId, jNetId);
			// 	inArc.weight = cut.get(key);
			// }
			auto jNetId = netArcs[inArc.decision].headId;
			auto key = getKey(qNetId, iNetId, jNetId);
			inArc.weight = cut.get(key);
			// auto jNetId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
			// inArc.weight = (inArc.decision != -1) ? cut.get(iNetId, qNetId, jNetId) : 0;
			node.state2 = inArc.weight + parentNode.state2;
			if (node.state2 < -0.5) nodesToRemove.push_back(nodeId);
		};
		for_each(layer.begin(), layer.end(), processNode);

		// if (!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove); // can use this to exit early.
	};

	// at this point, nodesToRemove contain ids of nodes from last layer, that are to be removed.

	size_t offset = tree.size()-1;
	for_each(tree.begin()+1, tree.begin()+offset, processLayer);

	// if entire tree need to be removed, exit from this tree.
	if (nodesToRemove.size() == tree[tree.size()-2].size()) return false;
	if (!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove);
	return true;
}

/******************************* END: RESTRICTED TREE *********************************************/



/******************************* BEGIN: RELAXED TREE **********************************************/


/**
 * Deletes the subtree of the given node recursively. Removes child nodes that have one incoming arc.
 */
void Inavap::RelaxedDD::topDownDelete(uint id) {
	/* start from the current node and recursively delete the children nodes until next node is
	 * terminal node or node with multiple incoming arcs. remove the last outgoing arcs of the last node.
	 */
	auto &node = nodes[id];
	auto arcsToDelete = node.outgoingArcs;

	for (auto outArcId: arcsToDelete) {
		// for each outArcId, find its head and apply topDownDelete() if head has single incoming arc.
		auto& outArc = arcs.at(outArcId);
		auto childId = outArc.head;
		auto &childNode = nodes[childId];
		deleteArc(node, outArc, childNode);
		if (childNode.incomingArcs.empty()) topDownDelete(childId); // orphan node
	}
	deleteNode(node);
}

/**
 * Deletes the node recursively until the root. Deletes the given node from the node's container.
 *
 * Function only deletes the ancestor nodes that doesn't have any children, called bereaved parent.
 */
void Inavap::RelaxedDD::bottomUpDelete(uint id) {

	// INVARIANT this node might contain multiple incoming parents, but should not contain any children.

	// for each incoming arc, remove arc and call soft delete on incoming node (iff has single outgoing arc).
	auto& node = nodes[id];
	auto incomingArcs = node.incomingArcs;
	for (const auto& arcId : incomingArcs){
		auto& arc = arcs[arcId];
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		//deleteArcById(arcId); // delete this incoming arc.
		if (parentNode.outgoingArcs.empty()){
			bottomUpDelete(parentNode.id);
		}
	}
	// delete this node from nodes.
	deleteNode(node);
}

void Inavap::RelaxedDD::deleteNode(LDDNode &node) {
	deletedNodeIds.emplace_back(node.id);
	nodes.erase(node.id);
}

[[always_inline]] void Inavap::RelaxedDD::deleteArc(LDDNode &tailNode, DDArc &arc, LDDNode &headNode) {
	// erase-remove idiom is replaced with erase(). below functions are available since C++20.
	erase(tailNode.outgoingArcs, arc.id);
	erase(headNode.incomingArcs, arc.id);
	arcs.erase(arc.id);
}

/**
 * Builds next layer from the given current layer. Returns the next layer containing the Ids of the nodes of next layer.
 * @param currentLayer
 * @param nextLayerSize
 * @param stateChangesNext
 * @return
 */
vector<uint> Inavap::RelaxedDD::buildNextLayer(const vector<uint> &currentLayer,
						uint& nextLayerSize, uint8_t stateChangesNext) {

	vector<uint> nextLayer;
	nextLayer.reserve(std::max(static_cast<uint>(4), nextLayerSize));
	nextLayerSize = 0;

	if (stateChangesNext) {
		/* next layer will undergo state change, so create only one node and do not set its state,
		 * its going to be changed in the next iteration anyway */

		LDDNode node {++lastInserted};

		for (auto id: currentLayer) {
			auto& parent = nodes[id];
			node.nodeLayer = parent.nodeLayer + 1;
			node.globalLayer = parent.globalLayer + 1;
			for (auto decision: parent.states) {
				// at least one arc should be matched to V-bar node. so remove all -1's in the matching.
				if (parent.states.size() > 1 && decision == -1) continue;
				auto nextArcId = ++lastInserted;
				DDArc newArc{nextArcId, id, node.id, decision};
				parent.outgoingArcs.push_back(nextArcId);
				node.incomingArcs.push_back(nextArcId);
				arcs.insert(make_pair(nextArcId, newArc));
			}
		}
		nextLayer.push_back(node.id);
		nodes.insert(make_pair(node.id, node));
	}
	else {

		/* The layers in the relaxed tree undergo state reduction. So all the states in the given layer are unique.
		 * Collect all the states in the (un-)ordered set and compare the membership of new state vector with the
		 * existing ones in the set. */
		// new strategy: do not do state reduction. build complete layer.
		for (auto id : currentLayer) {
			auto& parent = nodes[id];
			auto statesCopy = parent.states;

			for (auto decision : parent.states) {
				auto newStates(statesCopy);
				if (decision != -1)
					newStates.erase(find(newStates.begin(), newStates.end(), decision));
				auto nextId = ++lastInserted;
				LDDNode child{nextId};
				DDArc arc {nextId, id, nextId, decision};
				child.states = newStates;
				parent.outgoingArcs.push_back(nextId);
				child.incomingArcs.push_back(nextId);
				child.nodeLayer = parent.nodeLayer+1;
				child.globalLayer = parent.globalLayer+1;
				arcs.insert(make_pair(nextId, arc));
				nodes.insert(make_pair(nextId, child));
				nextLayer.push_back(nextId);
				nextLayerSize += child.states.size();
			}
		}
		// vector<LDDNode> nodesVector;
		// // vector<tuple<vector<int16_t>, int,int>> allStatesVector;
		// // unordered_set<tuple<set<int16_t>, int, int>, tuple_hash, tuple_equal> allStates;
		// // int j = 0;
		//
		// for (const auto id: currentLayer) {
		// 	auto &node = nodes[id];
		// 	auto statesCopy = node.states;
		//
		// 	for (auto decision: node.states) {
		// 		auto newStates(statesCopy);
		// 		if (decision != -1) newStates.erase(find(newStates.begin(), newStates.end(), decision));
		// 		auto nextId = ++lastInserted;
		//
		// 		// if newStates already in allStates, update exising node in nodesVector.
		// 		// set<int16_t> statesSet{newStates.begin(), newStates.end()};
		// 		// auto [it, isInserted] = allStates.insert({statesSet, 0, j});
		// 		// if (isInserted) { // create new Node
		// 			LDDNode newNode{nextId};
		// 			DDArc newArc{nextId, id, nextId, decision};
		// 			newNode.nodeLayer = node.nodeLayer+1;
		// 			newNode.globalLayer = node.globalLayer+1;
		// 			newNode.states = newStates;
		// 			newNode.incomingArcs.push_back(nextId);
		// 			node.outgoingArcs.push_back(nextId);
		// 			arcs.insert(make_pair(nextId, newArc));
		// 			nodes.insert(make_pair(node.id, newNode));
		// 			nodesVector.push_back(newNode);
		// 			// j++;
		// 		// } else { // state already exists, update existing node in nodesVector.
		// 		// 	auto [tempState, state2, pos] = *(it);
		// 		// 	auto &prevNode = nodesVector[pos];
		// 		// 	DDArc newArc{nextId, id, prevNode.id, decision};
		// 		// 	prevNode.incomingArcs.push_back(nextId);
		// 		// 	node.outgoingArcs.push_back(nextId);
		// 		// 	arcs.insert(make_pair(nextId, newArc));
		// 		// }
		// 	}
		// }
		//
		// // populate nextLayer. LATER. move nodes from vector.
		// for (auto& node : nodesVector) {
		// 	nextLayer.push_back(node.id);
		// 	nextLayerSize += node.states.size();
		// 	// nodes.insert(make_pair(node.id, node));
		// }
	}
	return nextLayer;
}

void Inavap::RelaxedDD::buildTree(Node node) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& stateUpdateMap = networkPtr->stateUpdateMap;

	// create root node and insert it.
	LDDNode root{node};
	startTree = node.globalLayer;
	// root.globalLayer = node.globalLayer;
	// // startTree = node.globalLayer;
	// root.nodeLayer = 0;
	// root.states = std::move(node.states);
	rootSolution = std::move(node.solutionVector);
	// insert root to the container.
	nodes.insert(make_pair(0, root));
	// setup root node and root layer.
	vector<uint> currentLayer = {0};
	tree.push_back(currentLayer);

	auto start = processingOrder.begin() + startTree;
	auto end = processingOrder.end();
	auto diff = std::distance(start, end);
	tree.reserve(diff+2); // root layer + terminal layer.

	uint nextLayerSize = 0;

	for (; start != end; ++start) {
		auto [a,b] = *start;
		// If V-bar changes in the next layer, states of nodes should be changed.
		if (stateUpdateMap.contains(a)) {
			const auto& newStates = stateUpdateMap.at(a);
			vector<int16_t> states;
			// remove this after changing in network.
			for (auto state: newStates) {
				states.push_back(static_cast<int16_t>(state));
			}
			updateStates(currentLayer, states);
		}

		/* at present, the current layer is copying the vector while inserting to the tree.
		 * Move current layer to the tree and get a pointer to the recently inserted layer
		 * of the tree. */

		vector<uint> nextLayer = buildNextLayer(currentLayer, nextLayerSize, networkPtr->hasStateChanged[a+1]);
		// reduce layer?
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}

	// build terminal node and terminal layer.
	vector<uint> terminalLayer;
	LDDNode terminalNode{++lastInserted};
	terminalId = terminalNode.id;

	/* the node layer, global layer, states, state2, and outgoing arcs are not relevant for the terminal node. */

	// create incoming arcs to terminal node. current layer now points to the last layer in the tree.
	for (auto id: currentLayer) {
		auto& parentNode = nodes[id];
		DDArc arc {++lastInserted, id, terminalId, 1};
		arc.weight = DOUBLE_MAX;
		terminalNode.incomingArcs.push_back(arc.id);
		parentNode.outgoingArcs.push_back(arc.id);
		arcs.insert(make_pair(arc.id, arc));
	}

	nodes.insert(make_pair(terminalId, terminalNode));
	terminalLayer.push_back(terminalId);
	tree.push_back(terminalLayer);
	// do not return cutset. might change later.
}

void Inavap::RelaxedDD::updateStates(const vector<uint>& currentLayer, const vector<int16_t>& nextLayerStates) {
	for (auto id: currentLayer) {
		auto& node = nodes[id];
		node.states = nextLayerStates;
	}
}

/**
 * Removes node and its associated subtree and parent(s) from the container.
 * @param id
 * @param isBatch
 */
void Inavap::RelaxedDD::removeNode(uint id, bool isBatch) {

	/* Deletes the incoming arcs from the parent, updates parent's outgoing arcs. Deletes node's outgoing arcs
	 * completely and updates the child's incoming arcs. Top-down delete is called for each child recursively until
	 * all the eligible children are deleted. Bottom-up delete is called recursively on each eligible ancestor.*/

	auto& node = nodes[id];

	/* According to the way that feasibility cut is applied to relaxed tree right now, we only remove nodes from the
     * last layer. This will not have any performance benefit if we call top-down delete on the terminal node multiple
     * times without deleting it. Thus, we can skip calling top-down delete function after removing the outgoing arc.
     * Uncomment the below block if new strategy is developed in applying the feasibility cut that involves deleting
     * nodes from the upper layers of the tree.
     */
	// auto outArcs = node.outgoingArcs;
	// for (auto childArcId : outArcs){
	// 	auto& childArc = arcs[childArcId];
	// 	auto& child = nodes[childArc.head];
	// 	deleteArc(node, childArc, child);
	//	// recursively delete sub-tree.
	// 	if (child.incomingArcs.empty()) topDownDelete(child.id); // orphan node
	// 	deleteNode(node); // comment this if uncommented the above line.
	// }

	// the last layer nodes have one outgoing arc to terminal, delete it.
	auto& outgoingArc = arcs[node.outgoingArcs.back()];
	auto& terminalNode = nodes[terminalId];
	deleteArc(node, outgoingArc, terminalNode);

	auto incomingArcs = node.incomingArcs;
	for (auto arcId: incomingArcs){
		auto& arc = arcs[arcId];
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		if (parentNode.outgoingArcs.empty()) bottomUpDelete(parentNode.id); // orphan node
	}
	// actual delete.
	deleteNode(node);

	if (!isBatch) updateTree();	// no batch deletion.
}

void Inavap::RelaxedDD::batchRemoveNodes(vector<uint> &nodes) {

	for (auto id: nodes) removeNode(id, false);
	updateTree();
}

/**
* Removes the deleted node Ids from the tree. Resets the deletedNodeIds variable after updating the tree. Optimized
* to call the function when removing nodes in batch.
*/
void Inavap::RelaxedDD::updateTree() {

	/* Since the node Ids are ordered by layers, the same code that is used in the restricted tree is used here.
	 * rewrite this function to if that invariant changes or the way the refinement applies to the relaxed tree changes.
	 */

	// sort by index (indirectly sorts by layer).
	sort(deletedNodeIds.begin(), deletedNodeIds.end());
	// uint count = deletedNodeIds.size();

	/* since both the deleted ids and layer is sorted, existence of an deleted Id in the given layer
	 * can be performed in O(1) time : compare with last element of layer.
	 *
	 * Invariant: current max in deleted Ids <= max Id in the layer.
	 */
	auto f = [this](vui& ids) {
		if (this->deletedNodeIds.empty()) return;
		//  vector doesn't support erasing with reverse iterator. create a mask for elements that are removed.
		vector<uint8_t> mask;
		mask.reserve(ids.size());

		for_each(ids.rbegin(), ids.rend(), [this,&mask](uint id) {
			// NOTE: the above invariant must be hold all time in order trick to work.
			if (!deletedNodeIds.empty() && id==this->deletedNodeIds.back()) {
				mask.push_back(1);
				this->deletedNodeIds.pop_back();
				// count--;
			}
			else mask.push_back(0); // edge case handled here.
		});

		std::reverse(mask.begin(), mask.end());
		size_t index = 0;
		ids.erase(std::remove_if(ids.begin(), ids.end(),
			[&mask,&index](uint id) { return static_cast<bool>(mask[index++]); }), ids.end());
	};
	for_each(tree.rbegin(), tree.rend(), f);
	deletedNodeIds.clear();
}
/**
 * Applies given optimality cut to the tree.
 * @param cut : optimality cut to apply on the tree
 * @return : upperbound of the relaxed DD.
 */
double Inavap::RelaxedDD::applyOptimalityCut(const Inavap::Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& rootNode = nodes[0];
	size_t i = 0;
	// compute justified RHS for root node.
	rootNode.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](const double val , const int decision) {
			if (decision == -1) {i++; return val;}
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			auto key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);

	size_t offset = tree.size()-1;
	size_t index = i; //startTree; // TODO set correct index.

	// for (size_t layer = 1; layer < tree.size()-1; layer++) {
	// 	const auto netArcId = processingOrder[i++].second;
	// 	const auto iNetId = netArcs[netArcId].tailId;
	// 	const auto qNetId = netArcs[netArcId].headId;
	//
	// 	auto processNode = [&cut, iNetId, qNetId, &netArcs, this](auto nodeId) {
	// 		auto& node = nodes.at(nodeId);
	// 		node.state2 = std::accumulate(node.incomingArcs.begin(), node.incomingArcs.end(), numeric_limits<double>::lowest(),
	// 			[this, &netArcs, &cut, &qNetId, &iNetId](const double val, const int inArcId) {
	// 				auto& inArc = arcs.at(inArcId);
	// 				inArc.weight = 0;
	// 				const auto& parentNode = nodes.at(inArc.tail);
	// 				if (inArc.decision == -1) return max(val, parentNode.state2); //return val;
	// 				// const auto& parentNode = nodes[inArc.tail];
	// 				auto jNetId = netArcs[inArc.decision].headId;
	// 				auto key = getKey(qNetId, iNetId, jNetId);
	// 				inArc.weight = cut.get(key);
	// 				return max(val, inArc.weight+parentNode.state2);
	// 		});
	// 	};
	// 	for_each(tree[layer].begin(), tree[layer].end(), processNode);
	// }

	/* Nested lambda functions to apply cut on each layer and each node in the selected layer respectively. The cut
	 * will be applied sequentially on each layer from the root to the terminal layer (both excluding). The functions
	 * will modify the node and its corresponding arcs while processing. */
	auto processLayer = [this, &processingOrder, &netArcs, &cut, &index](const auto& layer) {
		const auto netArcId = processingOrder[index++].second;
		const auto iNetId = netArcs[netArcId].tailId;
		const auto qNetId = netArcs[netArcId].headId;

		/* process each node sequentially from left to right, can use different execution policy (still sequential) */
		auto processNode =  [&cut, iNetId, qNetId, &netArcs, this](auto nodeId) {
			auto& node = nodes[nodeId];
			// double newState = std::numeric_limits<double>::lowest();
			// for (auto inArcId : node.incomingArcs) { // change to lambda function later.
			// 	auto& inArc = arcs[inArcId];
			// 	const auto& parentNode = nodes[inArc.tail];
			// 	auto jNetId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
			// 	inArc.weight = (inArc.decision != -1) ? cut.get(iNetId, qNetId, jNetId) : 0;
			// 	if (newState <= (inArc.weight + parentNode.state2)) newState = inArc.weight + parentNode.state2;
			// }
			// node.state2 = newState;

			/* This might be bit tricky to understand role of std::accumulate here. Iterate over the
			 * incoming arcs of the node and compute their weights and return maximum value. Instead
			 * of returning the accumulated value, return the max of computed arc weight and 'val'.
			 * If the decision on the arc is -1, return early.
			 */

			node.state2 = std::accumulate(node.incomingArcs.begin(), node.incomingArcs.end(), numeric_limits<double>::lowest(),
				[this, &netArcs, &cut, &qNetId, &iNetId](const double val, const int inArcId) {
					auto& inArc = arcs[inArcId];
					inArc.weight = 0;
					const auto& parentNode = nodes[inArc.tail];
					if (inArc.decision == -1) return max(val, parentNode.state2); //return val;
					// const auto& parentNode = nodes[inArc.tail];
					auto jNetId = netArcs[inArc.decision].headId;
					auto key = getKey(qNetId, iNetId, jNetId);
					inArc.weight = cut.get(key);
					return max(val, inArc.weight+parentNode.state2);
			});
		};
		for_each(layer.begin(), layer.end(), processNode); // call inner lambda function to process each node.
	};
	// for_each(tree.begin()+1, tree.begin()+offset, processLayer); // call outer lambda function to process layer.


	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		const auto netArcId = processingOrder[index++].second;
		const auto iNetId = netArcs[netArcId].tailId;
		const auto qNetId = netArcs[netArcId].headId;

		for (auto nodeId : tree[layer]) {
			auto& node = nodes[nodeId];
			// instead accumulate, use loop
			// get maximum state from all incoming arcs
			double maxState = DOUBLE_MIN;
			for (auto arcId : node.incomingArcs) {
				auto& inArc = arcs[arcId];
				inArc.weight = 0;
				const auto& parent = nodes[inArc.tail];
				if (inArc.decision != -1) {
					// get coeff from cut
					auto jNetId = netArcs[inArc.decision].headId;
					auto key = getKey(qNetId, iNetId, jNetId);
					inArc.weight = cut.get(key);
				}
				if ((inArc.weight + parent.state2) > maxState) {
					maxState = inArc.weight+parent.state2;
				}
			}
			node.state2 = maxState;
			// node.state2 = std::accumulate(node.incomingArcs.begin(), node.incomingArcs.end(), DOUBLE_MIN,
			// 	[this, &netArcs, &cut, &qNetId, &iNetId](const double val, const int inArcId) {
			// 		auto& inArc = arcs[inArcId];
			// 		inArc.weight = 0;
			// 		const auto& parentNode = nodes[inArc.tail];
			// 		if (inArc.decision == -1) return max(val, parentNode.state2); //return val;
			// 		// const auto& parentNode = nodes[inArc.tail];
			// 		auto jNetId = netArcs[inArc.decision].headId;
			// 		auto key = getKey(qNetId, iNetId, jNetId);
			// 		inArc.weight = cut.get(key);
			// 		return max(val, inArc.weight+parentNode.state2);
			// });
		}
	}

	// terminal layer.
	auto& terminalNode = nodes[terminalId];
	double terminalState = DOUBLE_MIN; // changing it to double::lowest() is acting weird, fix this later.

	/* From the incoming arcs of terminal node, update the arc with new state. Find the arc with maximum state. */
	for_each(terminalNode.incomingArcs.begin(), terminalNode.incomingArcs.end(),[this, &terminalState](auto arcId) {
		auto& arc = arcs[arcId];
		const auto& parentNode = nodes[arc.tail];
		arc.weight = std::min(arc.weight, parentNode.state2);
		terminalState = std::max(terminalState, arc.weight);
	});
	terminalNode.state2 = terminalState;
	return terminalState;
}

uint8_t Inavap::RelaxedDD::applyFeasibilityCut(const Inavap::Cut &cut) {
	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& rootNode = nodes[0];
	size_t i = 0;
	// compute justified RHS for root node.
	rootNode.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](double val , int decision) { // pass reference to val?
			if (decision == -1) {i++; return val; } // might not
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			uint64_t key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);

	size_t offset = tree.size()-1;
	size_t index = i;// startTree; // index of element in processing order vector. TODO verify this.

	// auto processLayer = [this, &processingOrder, &netArcs, &cut, &index](const auto& layer) {
	// 	const auto netArcId = processingOrder[index++].second;
	// 	const auto iNetId = netArcs[netArcId].tailId;
	// 	const auto qNetId = netArcs[netArcId].headId;
	//
	// 	auto processNode =  [&cut, iNetId, qNetId, &netArcs, this](auto nodeId) {
	// 		auto& node = nodes[nodeId];
	//
	// 		/* This might be bit tricky to understand role of std::accumulate here. Iterate over the
	// 		 * incoming arcs of the node and compute their weights and return maximum value. Instead
	// 		 * of returning the accumulated value, return the max of computed arc weight and 'val'.
	// 		 */
	//
	// 		node.state2 = std::accumulate(node.incomingArcs.begin(), node.incomingArcs.end(), numeric_limits<double>::lowest(),
	// 			[this, &netArcs, &cut, iNetId, qNetId](auto val, int inArcId) {
	// 				auto& inArc = arcs[inArcId];
	// 				inArc.weight = 0;
	// 				const auto &parentNode = nodes[inArc.tail];
	// 				if (inArc.decision == -1) return max(val, parentNode.state2);
	// 				// const auto& parentNode = nodes[inArc.tail];
	// 				auto jNetId = netArcs[inArc.decision].headId;
	// 				auto key = getKey(qNetId, iNetId, jNetId);
	// 				inArc.weight = cut.get(key);
	// 				return max(val, inArc.weight+parentNode.state2);
	// 			}
	// 		);
	// 		// auto& inArc = arcs[node.incomingArc];
	// 		// const auto& parentNode = nodes[inArc.tail];
	// 		// auto jNetId = (inArc.decision != -1) ? netArcs[inArc.decision].headId : 0;
	// 		// inArc.weight = (inArc.decision != -1) ? cut.get(iNetId, qNetId, jNetId) : 0;
	// 		// node.state2 = inArc.weight + parentNode.state2;
	// 	};
	// 	for_each(layer.begin(), layer.end(), processNode);
	// };
	//
	// for_each(tree.begin()+1, tree.begin()+offset, processLayer);

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		const auto netArcId = processingOrder[index++].second;
		const auto iNetId = netArcs[netArcId].tailId;
		const auto qNetId = netArcs[netArcId].headId;

		for (auto nodeId : tree[layer]) {
			auto& node = nodes[nodeId];
			double maxState = DOUBLE_MIN;
			for (auto arcId : node.incomingArcs) {
				auto& inArc = arcs[arcId];
				inArc.weight = 0;
				const auto& parent = nodes[inArc.tail];
				if (inArc.decision != -1) {
					// get coeff from cut
					auto jNetId = netArcs[inArc.decision].headId;
					auto key = getKey(qNetId, iNetId, jNetId);
					inArc.weight = cut.get(key);
				}
				if ((inArc.weight + parent.state2) > maxState) {
					maxState = inArc.weight+parent.state2;
				}
			}
			node.state2 = maxState;
			// node.state2 = std::accumulate(node.incomingArcs.begin(), node.incomingArcs.end(), numeric_limits<double>::lowest(),
			// 	[this, &netArcs, &cut, iNetId, qNetId](double val, int inArcId) {
			// 		auto& inArc = arcs[inArcId];
			// 		inArc.weight = 0;
			// 		const auto &parentNode = nodes[inArc.tail];
			// 		if (inArc.decision == -1) return max(val, parentNode.state2);
			// 		// const auto& parentNode = nodes[inArc.tail];
			// 		auto jNetId = netArcs[inArc.decision].headId;
			// 		auto key = getKey(qNetId, iNetId, jNetId);
			// 		inArc.weight = cut.get(key);
			// 		return max(val, inArc.weight+parentNode.state2);
			// 	}
			// );
		}
	}

	// nodes to remove.
	vui nodesToRemove;
	auto& lastLayer = tree[tree.size()-2];
	// for_each(lastLayer.begin(), lastLayer.end(), [this, &nodesToRemove](uint id) {
	// 	if (nodes[id].state2 < 0) nodesToRemove.push_back(id);
	// });
	for (auto nodeId : lastLayer) {
		if (nodes[nodeId].state2 < 0) nodesToRemove.push_back(nodeId);
	}

	if (nodesToRemove.size() == tree[tree.size()-2].size()) return false;
	if (!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove);
	return true;
}

/******************************************************************************************/

optional<vector<Inavap::Node>> Inavap::RestrictedDDNew::compile(Inavap::Node node) {
	const auto& processingOrder = networkPtr->processingOrder;
	const auto& stateUpdateMap = networkPtr->stateUpdateMap;
	const auto& netArcs = networkPtr->networkArcs;

	startTree = node.globalLayer;
	rootSolution = node.solutionVector;
	// create root node.
	RDDNode root{node}; // all the other attributes should be updated in constructor.
	nodes.insert(make_pair(0, root));

	Layer currentLayer;
	currentLayer.push_back(0);
	tree.push_back(currentLayer);

	auto start = processingOrder.begin() + startTree;
	auto end = processingOrder.end();
	uint index = 0;
	uint nextLayerSize = 0;
	uint16_t exactLayer = 0;
	uint8_t isExact = true;

	for (; start != end; ++start) {
		auto [a,b] = *start;
		if (stateUpdateMap.contains(a)) {
			const auto& newStates = stateUpdateMap.at(a);
			States tempStates;
			for (const auto s: newStates) tempStates.push_back(static_cast<int16_t>(s));
			// updateStates(currentLayer, tempStates);
			for_each(currentLayer.begin(), currentLayer.end(),
				[this, &tempStates](const auto id) {
						auto& node = nodes[id];
						node.states = tempStates;
			});
		}

		uint8_t stateChangesNext = networkPtr->hasStateChanged[a+1];
		Layer nextLayer = buildNextLayer(currentLayer, nextLayerSize, isExact,
				stateChangesNext);
		if (isExact) exactLayer++;
		++index;
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}

	terminalId = ++lastInserted;
	// need to update only id and state fields for terminal node.
	RDDNode terminal{terminalId};

	/* currentLayer should point to last layer of tree. Create terminal arc with decision '1'
	 * for every node in last layer.
	 */
	for (const auto id : currentLayer) {
		uint arcId = ++lastInserted;
		DDArc arc {arcId, id, terminalId, 1};
		arc.weight = std::numeric_limits<double>::max();
		nodes[id].outgoingArcs.push_back(arcId);
		terminalInArcs.push_back(arcId);
		arcs.insert(make_pair(arcId, arc));
	}

	Layer terminalLayer;
	terminalLayer.push_back(terminalId);
	tree.push_back(terminalLayer);
	nodes.insert(make_pair(terminalId, terminal));
	if (!isExact) status = 1;

	if (isExact) return nullopt;
	return getExactCutSet(tree[exactLayer]);
}

Inavap::Layer Inavap::RestrictedDDNew::buildRestrictedLayer(const Layer& currentLayer, uint8_t stateChangesNext) {

	Layer nextLayer;
	nextLayer.reserve(max_width);

	// if (stateChangesNext) {
	// 	for (const auto id: currentLayer) {
	// 		auto& node = nodes[id];
	// 		int16_t decision = -1;
	// 		for (auto start = node.states.rbegin(); start != node.states.rend(); ++start) {
	// 			auto state = *start;
	// 			if (state > decision ) decision = state;
	// 		}
	// 		// create child.
	// 		auto index = ++lastInserted;
	// 		node.outgoingArcs.push_back(index);
	// 		DDArc arc{index, id, index, decision};
	// 		RDDNode child{index};
	// 		child.incomingArc = index;
	// 		child.states = node.states;
	// 		// child.nodeLayer = node.nodeLayer+1;
	// 		// child.globalLayer = node.globalLayer + 1;
	// 		if (decision != -1) child.states.erase(find(child.states.begin(),
	// 			child.states.end(), decision));
	// 		nodes.insert(make_pair(index, child));
	// 		arcs.insert(make_pair(index, arc));
	// 		nextLayer.push_back(index);
	// 	}
	// 	return nextLayer;
	// }

	// auto layerNumber = nodes[currentLayer.front()].globalLayer;
	// if (networkPtr->hasStateChanged[layerNumber+1]) {
	// 	for (auto id : currentLayer) {
	// 		auto& node = nodes[id];
	// 		int16_t decision = -1;
	// 		for (auto start = node.states.rbegin(); )
	// 	}
	// }
	for (const auto id: currentLayer) {
		auto& node = nodes[id];
		// for each node, create only one child node.
		// auto decision = node.states.back(); // assuming states are sorted in asc order.
		auto decision = *std::max_element(node.states.begin(), node.states.end());
		auto index = ++lastInserted;
		node.outgoingArcs.push_back(index);
		DDArc arc{index, id, index, decision};
		RDDNode child{index};
		child.incomingArc = index;
		child.states = node.states;
		child.nodeLayer = node.nodeLayer+1;
		child.globalLayer = node.globalLayer + 1;
		if (decision != -1) child.states.erase(find(child.states.begin(),
			child.states.end(), decision));
		nodes.insert(make_pair(index, child));
		arcs.insert(make_pair(index, arc));
		nextLayer.push_back(index);
	}
	return nextLayer;
}

Inavap::Layer Inavap::RestrictedDDNew::buildNextLayer(const Layer& currentLayer,
	uint& nextLayerSize, uint8_t &isExact, uint8_t stateChangesnext) {

	if (isExact) {
		Layer nextLayer;
		nextLayer.reserve(nextLayerSize);
		nextLayerSize = 0;
		size_t count = 0;

		for (const auto id: currentLayer) {
			auto& node = nodes[id];
			const auto states = node.states;

			for (auto start = node.states.rbegin(); start != node.states.rend(); ++start) {
				auto decision = *start;

				if (count >= max_width) {isExact = false; return nextLayer;}

				auto index = ++lastInserted;
				node.outgoingArcs.push_back(index);
				RDDNode child{index};
				DDArc arc{index, id, index, decision};
				child.states = states;
				if (decision != -1) child.states.erase(find(child.states.begin(),
					child.states.end(), decision));
				child.incomingArc = index;
				child.globalLayer = node.globalLayer+1;
				child.nodeLayer = node.nodeLayer+1;
				nodes.insert(make_pair(index, child));
				arcs.insert(make_pair(index, arc));
				count++;
				nextLayer.push_back(index);
			}
		}
		return nextLayer;
	}
	// tree is not exact anymore
	return buildRestrictedLayer(currentLayer,stateChangesnext);
}

Inavap::Path Inavap::RestrictedDDNew::getPathForNode(uint id) const {

	Path solution;
	/* Since every node has single parent, recursively traverse until root */
	const RDDNode *current = &nodes.at(id);
	while (current->nodeLayer) { // current is not the root
		const auto& incomingArc = arcs.at(current->incomingArc);
		solution.push_back(incomingArc.decision);
		current = &nodes.at(incomingArc.tail);
	}
	// add root solution
	Path finalSolution{rootSolution};
	finalSolution.insert(finalSolution.end(), solution.rbegin(), solution.rend());
	return finalSolution;
}


vector<Inavap::Node> Inavap::RestrictedDDNew::getExactCutSet(const Layer& layer) const{
	vector<Inavap::Node> cutsetNodes;
	// for each node in the layer, get solution path and create Node object.
	for_each(layer.begin(), layer.end(), [&cutsetNodes,this](auto id) {
		const auto& node = nodes.at(id);
		cutsetNodes.emplace_back(node.states, getPathForNode(id), DOUBLE_MIN,
			DOUBLE_MIN, node.globalLayer);
	});
	return cutsetNodes;
}

Inavap::Path Inavap::RestrictedDDNew::getMaxPath() const {

	// ASAP update default lower bounds and upper bounds to something new
	uint maxId = 0; // .
	double maxWeight = DOUBLE_MIN;

	for (const auto arcId : terminalInArcs) {
		const auto& incomingArc = arcs.at(arcId);
		if (incomingArc.weight > maxWeight) {
			maxWeight = incomingArc.weight;
			maxId = incomingArc.tail;
		}
	}

	return getPathForNode(maxId);
}

void Inavap::RestrictedDDNew::removeNode(uint id, bool isBatch) {

	// for now only remove the id from terminalInArcs vector.
	auto& node = nodes[id];
	auto outArcId = node.outgoingArcs[0]; // should only have one out arc.
	auto res = std::find(terminalInArcs.begin(), terminalInArcs.end(), outArcId);
	if (res != terminalInArcs.end())
	terminalInArcs.erase(res);
	else cout << "Attempting to remove arc not present in the vector." << endl;
	node.outgoingArcs.clear();
	deletedNodeIds.push_back(id);
	nTerminalArcsRemoved++;

}

void Inavap::RestrictedDDNew::updateTree() {

	// for now only remove node ids from last layer of tree.
	auto& currentLayer = tree[tree.size()-2];
	for (auto id : deletedNodeIds) {
		currentLayer.erase(find(currentLayer.begin(), currentLayer.end(), id));
	}
	deletedNodeIds.clear();
}

void Inavap::RestrictedDDNew::batchRemoveNodes(const Layer& nodeIds) {

	for (auto id: nodeIds) {
		removeNode(id);
	}
	updateTree();
}

uint8_t Inavap::RestrictedDDNew::applyFeasibilityCut(const Inavap::Cut &cut) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& root = nodes[0];

	size_t i = 0;
	double justifiedRHS = cut.getRHS();

	root.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](double val , int decision) { // pass reference to val?
			if (decision == -1) {i++; return val; } // might not
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			uint64_t key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);

	// for (auto decision : rootSolution) { // replace with std::accumulate.
	// 	if (decision == -1) { i++; continue; }
	// 	auto netArcId = processingOrder[i++].second;
	// 	uint iNetId = netArcs[netArcId].tailId;
	// 	uint qNetId = netArcs[netArcId].headId;
	// 	uint jNetId = netArcs[decision].headId;
	// 	auto key = getKey(qNetId, iNetId, jNetId);
	// 	justifiedRHS += cut.get(key);
	// }

	// root.state2 = justifiedRHS;

	Layer nodesToRemove;

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArcId = processingOrder[i++].second;
		auto iNetId = netArcs[netArcId].tailId;
		auto qNetId = netArcs[netArcId].headId;

		nodesToRemove.clear();

		for (auto nodeId : tree[layer]) {
			auto& node = nodes[nodeId];
			auto& inArc = arcs[node.incomingArc];
			const auto& parent = nodes[inArc.tail];

			if (inArc.decision != -1) {
				auto jNetId = netArcs[inArc.decision].headId;
				auto key = getKey(qNetId, iNetId, jNetId);
				inArc.weight = cut.get(key);
				node.state2 = parent.state2 + inArc.weight;
			}
			else {
				inArc.weight = 0.0;
				node.state2 = parent.state2;
			}
			if (node.state2 < -0.5) nodesToRemove.push_back(nodeId);
		}
		// somehow this lambda function is causing the bug. inspect this.
		// for_each(std::execution::seq, tree[layer].begin(), tree[layer].end(),
		// 	[this,&netArcs, &cut, iNetId, qNetId, &nodesToRemove](const auto id) {
		// 			auto& node = nodes[id];
		// 			auto& inArc = arcs[node.incomingArc];
		// 			const auto& parent = nodes[inArc.tail];
		//
		// 			if (inArc.decision == -1) {node.state2 = parent.state2; return;}
		// 			auto jNetId = netArcs[inArc.decision].headId;
		// 			auto key = getKey(qNetId, iNetId, jNetId);
		// 			inArc.weight = cut.get(key);
		// 			node.state2 = inArc.weight + parent.state2;
		// 			if (node.state2 < -0.5) nodesToRemove.push_back(id);
		// 		}
		// );
		// for_each(std::execution::seq, tree[layer].begin(), tree[layer].end(), processNode);
	}

	// only terminal arcs are removed during node removal.
	if (terminalInArcs.empty()) return false; // later return to previous version.
	if (!nodesToRemove.empty()) batchRemoveNodes(nodesToRemove);
	if (terminalInArcs.empty()) return false;
	return true;
}

double Inavap::RestrictedDDNew::applyOptimalityCut(const Inavap::Cut &cut) {
	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	auto& root = nodes[0];
	double justifiedRHS = cut.getRHS();
	size_t i = 0;
	root.state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
[&processingOrder, &netArcs, &i, &cut](double val , int decision) { // pass reference to val?
			if (decision == -1) {i++; return val; } // might not
			auto netArcId = processingOrder[i++].second;
			auto iNetId = netArcs[netArcId].tailId;
			auto qNetId = netArcs[netArcId].headId;
			auto jNetId = netArcs[decision].headId;
			uint64_t key = getKey(qNetId, iNetId, jNetId);
			return val + cut.get(key);
		}
	);
	// for (auto decision: rootSolution) {
	// 	if (decision == -1) {i++; continue;}
	// 	auto netArcId = processingOrder[i++].second;
	// 	auto iNetId = netArcs[netArcId].tailId;
	// 	auto qNetId = netArcs[netArcId].headId;
	// 	uint jNetId = netArcs[decision].headId;
	// 	auto key = getKey(qNetId, iNetId, jNetId);
	// 	justifiedRHS += cut.get(key);
	// }
	// root.state2 = justifiedRHS;

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArcId = processingOrder[i++].second;
		auto iNetId = netArcs[netArcId].tailId;
		auto qNetId = netArcs[netArcId].headId;

		for_each(std::execution::seq, tree[layer].begin(), tree[layer].end(),
			[this, &netArcs, &cut, iNetId, qNetId](const auto id) {
					auto& node = nodes[id];
					auto& inArc = arcs[node.incomingArc];
					const auto& parent = nodes[inArc.tail];
					if (inArc.decision == -1) {node.state2 = parent.state2; return;}
					auto jNetId = netArcs[inArc.decision].headId;
					auto key = getKey(qNetId, iNetId, jNetId);
					inArc.weight = cut.get(key);
					node.state2 = inArc.weight + parent.state2;
			}
		);

		// for (auto nodeId : tree[layer]) {
		// 	auto& node = nodes[nodeId];
		// 	auto inArc = arcs[node.incomingArc];
		// 	const auto& parent = nodes[inArc.tail];
		// 	// if (find(parent.states.begin(), parent.states.end(), inArc.decision)
		// 	// 		!= parent.states.end()) {
		// 		if (inArc.decision != -1) {
		// 			auto jNetId = netArcs[inArc.decision].headId;
		// 			auto key = getKey(qNetId, iNetId, jNetId);
		// 			inArc.weight = cut.get(key);
		// 			node.state2 = parent.state2 + inArc.weight;
		// 		}
		// 		else {
		// 			inArc.weight = 0.0;
		// 			node.state2 = parent.state2;
		// 		}
		// 	// }
		// }
	}

	// terminal arcs.
	auto& terminal = nodes[terminalId];
	double terminalState =  DOUBLE_MIN;
	for_each(std::execution::seq, terminalInArcs.begin(), terminalInArcs.end(),
		[this, &terminalState](const auto id) {
				auto& arc = arcs[id];
				const auto& parent = nodes[arc.tail];
				arc.weight = min(arc.weight, parent.state2);
				terminalState = max(terminalState, arc.weight);
		}
	);
	terminal.state2 = terminalState;
	return terminalState;
}

/* ******************************************************************************************************************** */

void Inavap::RelaxedDDNew::buildTree(Node node) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& stateUpdateMap = networkPtr->stateUpdateMap;

	startTree = node.globalLayer;
	rootSolution = std::move(node.solutionVector);

	// create root node and insert it.
	LDDNode root{node};
	nodes.insert(make_pair(0, root));
	// setup root node and root layer.
	// vector<uint> currentLayer = {0};
	tree.push_back({0});

	auto start = processingOrder.begin() + startTree;
	auto end = processingOrder.end();
	auto diff = std::distance(start, end);
	tree.reserve(diff+2); // root layer + terminal layer.

	uint nextLayerSize = 0;
	uint index = 0;

	for (index = 0; start != end; ++start, ++index) {
		auto [a,b] = *start;
		// If V-bar changes in the next layer, states of nodes should be changed.
		if (stateUpdateMap.contains(a)) {
			const auto& newStates = stateUpdateMap.at(a);
			vector<int16_t> states;
			// remove this after changing in network.
			for (auto state: newStates) {
				states.push_back(static_cast<int16_t>(state));
			}
			for_each(tree[index].cbegin(), tree[index].cend(),
				[this, &states](const auto id) {
						auto& node = nodes[id];
						node.states = states;
			});
			// states of current layer nodes is updated, nextLayerSize is not accurate.
			nextLayerSize = tree[index].size() * states.size();
		}
		uint8_t stateChangesNext = networkPtr->hasStateChanged[a+1];
		buildNextLayer(index, nextLayerSize, stateChangesNext);
	}

	// build terminal node and terminal layer.
	vector<uint> terminalLayer;
	LDDNode terminalNode{++lastInserted};
	terminalId = terminalNode.id;

	terminalNode.state2 = DOUBLE_MIN;
	/* the node layer, global layer, states, state2, and outgoing arcs are not relevant for the terminal node. */

	// create incoming arcs to terminal node. current layer now points to the last layer in the tree.
	for (auto id: tree[index]) { // is it index or index-1?
		auto& parentNode = nodes[id];
		DDArc arc {++lastInserted, id, terminalId, 1};
		arc.weight = DOUBLE_MAX;
		terminalNode.incomingArcs.push_back(arc.id);
		parentNode.outgoingArcs.push_back(arc.id);
		arcs.insert(make_pair(arc.id, arc));
	}

	nodes.insert(make_pair(terminalId, terminalNode));
	terminalLayer.push_back(terminalId);
	tree.push_back(terminalLayer);
}


void Inavap::RelaxedDDNew::buildNextLayer(uint current, uint &nextLayerSize, uint8_t stateChangesNext) {

	const auto& currentLayer = tree[current];

	/* If the next layer size (@param nextLayerSize) is greater than MAX WIDTH, then next layer undergo reduction
	 * to single node. So, create one node and make all the outgoing arcs from the nodes in the current layer point
	 * to the newly created node in the next layer.
	 *
	 * Otherwise, create nodes as usual (new node for every outgoing arc from current layer).
	 */
	// if (nextLayerSize >= RELAXED_MAX_WIDTH );
	if (nextLayerSize >= RELAXED_MAX_WIDTH && (nodes[currentLayer.front()].globalLayer < networkPtr->totalLayers-5)) {
		status = NON_EXACT; // set DD status flag to non-exact.

		LDDNode newNode {++lastInserted};
		newNode.incomingArcs.reserve(nextLayerSize); // reserve space for incoming arcs.

		set<int16_t> allStates;
		std::for_each(currentLayer.cbegin(), currentLayer.cend(), [&](uint id) {
			auto& node = nodes[id];
			const auto& states = node.states;
			allStates.insert(states.begin(), states.end()); // union of all states.
			// create new arc for each state.
			std::for_each(states.begin(), states.end(), [&](int16_t state) {
				if (true || !(stateChangesNext && node.states.size() > 1 && state == -1)) { // likely branch taken.
					DDArc newArc {++lastInserted, id, newNode.id, state};
					node.outgoingArcs.push_back(newArc.id);
					newNode.incomingArcs.push_back(newArc.id);
					arcs.insert(make_pair(newArc.id, newArc));
				}
			});
		});

		nextLayerSize = allStates.size(); // update nextLayerSize.

		const auto& someParent = nodes[currentLayer[0]];
		newNode.nodeLayer = someParent.nodeLayer + 1;
		newNode.globalLayer = someParent.globalLayer + 1;
		newNode.states = std::vector(allStates.begin(), allStates.end());
		// newNode.outgoingArcs.reserve(newNode.states.size()); // reserve space for outgoing arcs.
		nodes.insert(make_pair(newNode.id, newNode));
		tree.push_back({newNode.id});
		return ;
	}

	vector<uint> nextLayer;
	nextLayer.reserve(std::max(nextLayerSize, static_cast<uint>(4)));
	nextLayerSize = 0;

	std::for_each(currentLayer.begin(), currentLayer.end(), [&](uint id) {
		auto& parent = nodes[id];
		const auto& states = parent.states;

		std::for_each(states.begin(), states.end(), [&](int16_t state) {
			if (true || !(stateChangesNext && states.size() > 1 && state == -1)) {
				auto newStates = states;
				if (state != -1) // remove selected decision from the states.
					newStates.erase(std::remove(newStates.begin(), newStates.end(), state),
										newStates.end());

				nextLayerSize += newStates.size();
				auto nextId = ++lastInserted;
				DDArc childArc{nextId, id, nextId, state};
				LDDNode childNode{nextId};
				childNode.incomingArcs.push_back(nextId);
				childNode.nodeLayer = parent.nodeLayer+1;
				childNode.globalLayer = parent.globalLayer+1;
				childNode.states = std::move(newStates);
				parent.outgoingArcs.push_back(nextId);
				arcs.insert(make_pair(nextId, childArc));
				nodes.insert(make_pair(nextId, childNode));
				nextLayer.push_back(nextId);
			}
		});
	});

	tree.push_back(std::move(nextLayer));

	return ;

	if (stateChangesNext) {
		/* next layer will undergo state change, so create only one node and do not set its state,
		 * its going to be changed in the next iteration anyway */

		LDDNode node {++lastInserted};

		for (auto id: currentLayer) {
			auto& parent = nodes[id];
			node.nodeLayer = parent.nodeLayer + 1;
			node.globalLayer = parent.globalLayer + 1;
			for (auto decision: parent.states) {
				// at least one arc should be matched to V-bar node. so remove all -1's in the matching.
				if (parent.states.size() > 1 && decision == -1) continue;
				auto nextArcId = ++lastInserted;
				DDArc newArc{nextArcId, id, node.id, decision};
				parent.outgoingArcs.push_back(nextArcId);
				node.incomingArcs.push_back(nextArcId);
				arcs.insert(make_pair(nextArcId, newArc));
			}
		}
		nextLayer.push_back(node.id);
		nodes.insert(make_pair(node.id, node));
	}
	else {

		/* The layers in the relaxed tree undergo state reduction. So all the states in the given layer are unique.
		 * Collect all the states in the (un-)ordered set and compare the membership of new state vector with the
		 * existing ones in the set. */
		// new strategy: do not do state reduction. build complete layer.
		for (auto id : currentLayer) {
			auto& parent = nodes[id];
			auto statesCopy = parent.states;

			for (auto decision : parent.states) {
				auto newStates(statesCopy);
				if (decision != -1)
					newStates.erase(find(newStates.begin(), newStates.end(), decision));
				auto nextId = ++lastInserted;
				LDDNode child{nextId};
				DDArc arc {nextId, id, nextId, decision};
				child.states = newStates;
				parent.outgoingArcs.push_back(nextId);
				child.incomingArcs.push_back(nextId);
				child.nodeLayer = parent.nodeLayer+1;
				child.globalLayer = parent.globalLayer+1;
				arcs.insert(make_pair(nextId, arc));
				nodes.insert(make_pair(nextId, child));
				nextLayer.push_back(nextId);
				nextLayerSize += child.states.size();
			}
		}
		// vector<LDDNode> nodesVector;
		// // vector<tuple<vector<int16_t>, int,int>> allStatesVector;
		// // unordered_set<tuple<set<int16_t>, int, int>, tuple_hash, tuple_equal> allStates;
		// // int j = 0;
		//
		// for (const auto id: currentLayer) {
		// 	auto &node = nodes[id];
		// 	auto statesCopy = node.states;
		//
		// 	for (auto decision: node.states) {
		// 		auto newStates(statesCopy);
		// 		if (decision != -1) newStates.erase(find(newStates.begin(), newStates.end(), decision));
		// 		auto nextId = ++lastInserted;
		//
		// 		// if newStates already in allStates, update exising node in nodesVector.
		// 		// set<int16_t> statesSet{newStates.begin(), newStates.end()};
		// 		// auto [it, isInserted] = allStates.insert({statesSet, 0, j});
		// 		// if (isInserted) { // create new Node
		// 			LDDNode newNode{nextId};
		// 			DDArc newArc{nextId, id, nextId, decision};
		// 			newNode.nodeLayer = node.nodeLayer+1;
		// 			newNode.globalLayer = node.globalLayer+1;
		// 			newNode.states = newStates;
		// 			newNode.incomingArcs.push_back(nextId);
		// 			node.outgoingArcs.push_back(nextId);
		// 			arcs.insert(make_pair(nextId, newArc));
		// 			nodes.insert(make_pair(node.id, newNode));
		// 			nodesVector.push_back(newNode);
		// 			// j++;
		// 		// } else { // state already exists, update existing node in nodesVector.
		// 		// 	auto [tempState, state2, pos] = *(it);
		// 		// 	auto &prevNode = nodesVector[pos];
		// 		// 	DDArc newArc{nextId, id, prevNode.id, decision};
		// 		// 	prevNode.incomingArcs.push_back(nextId);
		// 		// 	node.outgoingArcs.push_back(nextId);
		// 		// 	arcs.insert(make_pair(nextId, newArc));
		// 		// }
		// 	}
		// }
		//
		// // populate nextLayer. LATER. move nodes from vector.
		// for (auto& node : nodesVector) {
		// 	nextLayer.push_back(node.id);
		// 	nextLayerSize += node.states.size();
		// 	// nodes.insert(make_pair(node.id, node));
		// }
	}
}

Inavap::Path Inavap::RelaxedDDNew::getPathForNode(uint id) const {

	// recursively traverse from current node till root.
	const LDDNode *current = &nodes.at(id);
	Path path;

	while (current->nodeLayer) {
		const LDDNode *potential_parent = &nodes.at(arcs.at(current->incomingArcs[0]).tail);
		//
		for (const auto arcId : current->incomingArcs) {
			const auto& arc = arcs.at(arcId);
			const LDDNode *parent = &nodes.at(arc.tail);
			if ((parent->state2 + arc.weight) == current->state2) {
				path.push_back(arc.decision);
				potential_parent = parent;
				break;
			}
		}
		current = potential_parent;
	}
	// prepend root solution.
	Path solution(rootSolution.begin(), rootSolution.end());
	solution.insert(solution.end(), path.rbegin(), path.rend());
	return solution;
}

/**
 * @return A maximum path.
 */
Inavap::Path Inavap::RelaxedDDNew::getSolution() const {

	uint maxId = 0; // .
	double maxWeight = DOUBLE_MIN;
	/* find a terminal arc with max weight and get solution for the incoming node of the arc.*/
	const auto& terminalNode = nodes.at(terminalId);
	std::for_each(terminalNode.incomingArcs.begin(), terminalNode.incomingArcs.end(),
		[&](auto id) {
			const auto& incomingArc = arcs.at(id);
			if (incomingArc.weight > maxWeight) {
				maxWeight = incomingArc.weight;
				maxId = incomingArc.tail;
			}
	});
	return getPathForNode(maxId);
}

uint8_t Inavap::RelaxedDDNew::applyFeasibilityCut(const Inavap::Cut &cut) {
	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for root node.
	size_t i = 0;
	nodes[0].state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
	[&processingOrder, &netArcs, &i, &cut](double val , int decision) {
                if (decision == -1) {i++; return val; } // might not
                auto netArcId = processingOrder[i++].second;
                auto iNetId = netArcs[netArcId].tailId;
                auto qNetId = netArcs[netArcId].headId;
                auto jNetId = netArcs[decision].headId;
                uint64_t key = getKey(qNetId, iNetId, jNetId);
                return val + cut.get(key);
			}
	);

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArcId = processingOrder[i++].second;
		auto iNetId = netArcs[netArcId].tailId;
		auto qNetId = netArcs[netArcId].headId;

		for (auto nodeId : tree[layer]) {
			auto& node = nodes[nodeId];
			double newState = DOUBLE_MIN;
			for (auto inArcId : node.incomingArcs) {
				auto& inArc = arcs[inArcId];
				const auto& parent = nodes[inArc.tail];
				if (inArc.decision != -1) {
                    auto jNetId = netArcs[inArc.decision].headId;
                    auto key = getKey(qNetId, iNetId, jNetId);
                    inArc.weight = cut.get(key);
					// no need to save arc weights at all. redundant store to inArc.weight;
                    newState = max(parent.state2 + inArc.weight, newState);
                }
				else newState = max(newState, parent.state2);
			}
			node.state2 = newState;
		}
	}

	vui nodesToRemove;
	size_t llayer = tree.size()-2;
	for (auto nodeId: tree[llayer]) {
		if (nodes[nodeId].state2 < -0.01)
			nodesToRemove.push_back(nodeId);
	}

	if (nodesToRemove.size() == tree[llayer].size()) return false;
	if (!nodesToRemove.empty())
		batchRemoveNodes(nodesToRemove);

	if (!isTreeExact()) {
		double maxState = DOUBLE_MIN;
		for (auto nodeId : tree[llayer]) {
			maxState = max(maxState, nodes[nodeId].state2);
		}

		vui arcsToRemoved;
		for (size_t layer = 1; layer < llayer; layer++) {
			if (tree[layer].size() ==1) {
				double maxGain = maxState - nodes[tree[layer][0]].state2;
				for (auto nodeId : tree[layer-1]) {
					auto& node = nodes[nodeId];
					for (auto outArcId : node.outgoingArcs) {
						auto& outArc = arcs[outArcId];
						const auto& parentNode = nodes[outArc.tail];
						if ((parentNode.state2 + outArc.weight + maxGain) <= 0.01) {
							arcsToRemoved.push_back(outArcId);
						}
					}
				}
			}
		}
		if (!arcsToRemoved.empty()) {
			batchRemoveArcs(arcsToRemoved);
		}
	}
	return true;
}

double Inavap::RelaxedDDNew::applyOptimalityCut(const Inavap::Cut &cut, double optimal, double upperbound) {

	const auto& processingOrder = networkPtr->processingOrder;
	const auto& netArcs = networkPtr->networkArcs;

	// compute justified RHS for root node.
	size_t i = 0;
	nodes[0].state2 = std::accumulate(rootSolution.begin(), rootSolution.end(), cut.getRHS(),
	[&processingOrder, &netArcs, &i, &cut](double val , int decision) {
                if (decision == -1) {i++; return val; } // might not
                auto netArcId = processingOrder[i++].second;
                auto iNetId = netArcs[netArcId].tailId;
                auto qNetId = netArcs[netArcId].headId;
                auto jNetId = netArcs[decision].headId;
                uint64_t key = getKey(qNetId, iNetId, jNetId);
                return val + cut.get(key);
			}
	);

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArcId = processingOrder[i++].second;
		auto iNetId = netArcs[netArcId].tailId;
		auto qNetId = netArcs[netArcId].headId;

		for (auto nodeId : tree[layer]) {
			auto& node = nodes[nodeId];
			double newState = DOUBLE_MIN;
			for (auto inArcId : node.incomingArcs) {
				auto& inArc = arcs[inArcId];
				const auto& parent = nodes[inArc.tail];
				if (inArc.decision != -1) {
                    auto jNetId = netArcs[inArc.decision].headId;
                    auto key = getKey(qNetId, iNetId, jNetId);
                    inArc.weight = cut.get(key);
					// no need to save arc weights at all. redundant store to inArc.weight;
                    newState = max(parent.state2 + inArc.weight, newState);
                }
				else newState = max(newState, parent.state2);
			}
			node.state2 = newState;
		}
	}

	double terminalState = DOUBLE_MIN;
	auto& terminalNode = nodes[terminalId];
	for (auto inArcId : terminalNode.incomingArcs) {
		auto& inArc = arcs[inArcId];
		const auto& parentNode = nodes[inArc.tail];
		// min of terminal arc weight, and parent state, and maximum of terminal arc weights
		inArc.weight = std::min(inArc.weight, parentNode.state2);
		terminalState = std::max(terminalState, inArc.weight);
	}
	terminalNode.state2 = terminalState;
	if (terminalState <= optimal) return terminalState;

	if (!isTreeExact()) {
		size_t llayer = tree.size() -2;
		double maxState = DOUBLE_MIN;
		for (auto nodeId : tree[llayer]) {
			maxState = max(maxState, nodes[nodeId].state2);
		}

		vui arcsToRemoved;
		for (size_t layer = 3; layer < llayer-1; layer++) {
			if (tree[layer].size() ==1) {
				double maxGain = maxState - nodes[tree[layer][0]].state2;
				for (auto nodeId : tree[layer-1]) {
					auto& node = nodes[nodeId];
					for (auto outArcId : node.outgoingArcs) {
						auto& outArc = arcs[outArcId];
						const auto& parentNode = nodes[outArc.tail];
						if ((parentNode.state2 + outArc.weight + maxGain) <= (optimal - 0.01)) {
							arcsToRemoved.push_back(outArcId);
						}
					}
				}
			}
		}
		if (!arcsToRemoved.empty()) {
			batchRemoveArcs(arcsToRemoved);
		}
	}
	return terminalState;
}

void Inavap::RelaxedDDNew::deleteNode(LDDNode &node) {
	deletedNodeIds.emplace_back(node.id);
	nodes.erase(node.id);
}

void Inavap::RelaxedDDNew::deleteArc(LDDNode &parent, DDArc &arc, LDDNode &child) {
	// erase-remove idiom is replaced with erase(). below functions are available since C++20.
	std::erase(parent.outgoingArcs, arc.id);
	std::erase(child.incomingArcs, arc.id);
	arcs.erase(arc.id);
}

void Inavap::RelaxedDDNew::removeNode(uint id, bool isBatch) {

	// this function is called on the nodes that are in last layer.
	/* Deletes the incoming arcs from the parent, updates parent's outgoing arcs. Deletes node's outgoing arcs
	 * completely and updates the child's incoming arcs. Top-down delete is called for each child recursively until
	 * all the eligible children are deleted. Bottom-up delete is called recursively on each eligible ancestor.*/

	auto& node = nodes[id];

	/* According to the way that feasibility cut is applied to relaxed tree right now, we only remove nodes from the
     * last layer. This will not have any performance benefit if we call top-down delete on the terminal node multiple
     * times without deleting it. Thus, we can skip calling top-down delete function after removing the outgoing arc.
     * Uncomment the below block if new strategy is developed in applying the feasibility cut that involves deleting
     * nodes from the upper layers of the tree.
     */
	// auto outArcs = node.outgoingArcs;
	// for (auto childArcId : outArcs){
	// 	auto& childArc = arcs[childArcId];
	// 	auto& child = nodes[childArc.head];
	// 	deleteArc(node, childArc, child);
	//	// recursively delete sub-tree.
	// 	if (child.incomingArcs.empty()) topDownDelete(child.id); // orphan node
	// 	deleteNode(node); // comment this if uncommented the above line.
	// }

	// the last layer nodes have one outgoing arc to terminal, delete it.
	auto& outgoingArc = arcs[node.outgoingArcs.back()];
	auto& terminalNode = nodes[terminalId];
	deleteArc(node, outgoingArc, terminalNode);

	auto incomingArcs = node.incomingArcs;
	for (auto arcId: incomingArcs){
		auto& arc = arcs[arcId];
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		if (parentNode.outgoingArcs.empty()) bottomUpDelete(parentNode.id); // orphan node
	}
	// actual delete.
	deleteNode(node);

	if (!isBatch) updateTree();	// no batch deletion.

}

void Inavap::RelaxedDDNew::bottomUpDelete(uint id) {
	// INVARIANT this node might contain multiple incoming parents, but should not contain any children.

	// for each incoming arc, remove arc and call soft delete on incoming node (iff has single outgoing arc).
	auto& node = nodes[id];
	auto incomingArcs = node.incomingArcs;
	for (const auto& arcId : incomingArcs){
		auto& arc = arcs[arcId];
		auto& parentNode = nodes[arc.tail];
		deleteArc(parentNode, arc, node);
		//deleteArcById(arcId); // delete this incoming arc.
		if (parentNode.outgoingArcs.empty()){
			bottomUpDelete(parentNode.id);
		}
	}
	// delete this node from nodes.
	deleteNode(node);
}

void Inavap::RelaxedDDNew::updateTree() {

	// sort by index (indirectly sorts by layer).
	sort(deletedNodeIds.begin(), deletedNodeIds.end());

	/* since both the deleted ids and layer is sorted, existence of an deleted Id in the given layer
	 * can be performed in O(1) time : compare with last element of layer.
	 *
	 * Invariant: current max in deleted Ids <= max Id in the layer.
	 */
	auto f = [this](vui& ids) {
		if (this->deletedNodeIds.empty()) return;
		//  vector doesn't support erasing with reverse iterator. create a mask for elements that are removed.
		vector<uint8_t> mask;
		mask.reserve(ids.size());

		for_each(ids.rbegin(), ids.rend(), [this,&mask](uint id) {
			// NOTE: the above invariant must be hold all time in order trick to work.
			if (!deletedNodeIds.empty() && id==this->deletedNodeIds.back()) {
				mask.push_back(1);
				this->deletedNodeIds.pop_back();
			}
			else mask.push_back(0); // edge case handled here.
		});

		std::reverse(mask.begin(), mask.end());
		size_t index = 0;
		ids.erase(std::remove_if(ids.begin(), ids.end(),
			[&mask,&index](uint id) { return static_cast<bool>(mask[index++]); }), ids.end());
	};
	for_each(tree.rbegin(), tree.rend(), f);
	deletedNodeIds.clear();
}

void Inavap::RelaxedDDNew::batchRemoveNodes(const vui &nodeIds) {
	for (auto nodeId : nodeIds) {
		removeNode(nodeId);
	}
	updateTree();
}

void Inavap::RelaxedDDNew::batchRemoveArcs(const vui &arcIds) {
	// calling deleteArcById on all arcsIds.
	for (auto arcId : arcIds) {
		auto& arc  = arcs[arcId];
		auto& tail = nodes[arc.tail];
		auto& head = nodes[arc.head];

		std::erase(head.incomingArcs, arcId);
		std::erase(tail.outgoingArcs, arcId);
		// head.incomingArcs.erase(std::find(head.incomingArcs.begin(), head.incomingArcs.end(), arcId));
		// tail.outgoingArcs.erase(std::find(tail.outgoingArcs.begin(), tail.outgoingArcs.end(), arcId));
		arcs.erase(arcId);

		// if (head.incomingArcs.empty()) removeNode(head.id);
		// if (tail.outgoingArcs.empty()) removeNode(tail.id);
	}
}

vector<Inavap::Node> Inavap::RelaxedDDNew::getCutset(double ub) {

	uint layer = 3;
	while (tree[layer].size() != 1) layer++; // advance to first non-exact layer.

	assertm(tree[layer].size() ==1, "This should not be cutset layer.");

	uint globalLayerNumber = nodes.at(tree[layer][0]).globalLayer;
	vector<Node> cutsetNodes;

	if (networkPtr->hasStateChanged[globalLayerNumber]) { // if state changes, get new states.
		const auto newStates = networkPtr->stateUpdateMap.at(globalLayerNumber);
		const vector<int16_t> statesVector (newStates.begin(), newStates.end());
		for (auto nodeId : tree[layer-1]) {
			const auto& node = nodes.at(nodeId);
			const auto partialSol = getPathForNode(nodeId);
			for (auto arcId : node.outgoingArcs) {
				auto decision = arcs.at(arcId).decision;
				auto solution = partialSol;
				solution.push_back(decision);
				cutsetNodes.emplace_back(statesVector, std::move(solution),DOUBLE_MIN, ub, globalLayerNumber);
			}
		}
	}
	else {
		for (auto nodeId : tree[layer-1]) {
			const auto& node = nodes.at(nodeId);
			const auto partialSol = getPathForNode(nodeId);
			for (auto arcId : node.outgoingArcs) {
				auto decision = arcs.at(arcId).decision;
				auto solution = partialSol;
				solution.push_back(decision);
				auto states = node.states;
				if (decision != -1) std::erase(states, decision);
				cutsetNodes.emplace_back(std::move(states), std::move(solution), DOUBLE_MIN, ub, globalLayerNumber);
			}
		}
	}
	return cutsetNodes;
}