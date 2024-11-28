//
// Created by nandgate on 6/3/24.
//


#include "DD.h"
#include <queue>
#include <limits>

/**
 * Compiles the decision diagram tree with given node as the root.
 */
optional<vector<Node_t>> DD::build(DDNode &node,int type_param) {
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
		isExact = buildNextLayer(currentLayer, nextLayer, networkPtr->hasStateChanged[a+1],type_param);
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
bool DD::buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, bool stateChangesNext,int type_param) {
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
			if(type_param == 1) {
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
			 }else {
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
			 			// 	auto [tempState, state2, pos] = *(it);
			 			// 	auto &prevNode = nodesVector[pos];
			 			// 	DDArc newArc{nextId, id, prevNode.id, decision};
			 			// 	prevNode.incomingArcs.push_back(nextId);
			 			// 	node.outgoingArcs.push_back(nextId);
			 			// 	arcs.insert(make_pair(nextId, newArc));
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
			else {
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
						 auto [it, isInserted] = allStates.insert({newStates, 0, j});
						 if (isInserted) { // create new Node
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
						 } else { // state already exists, update existing node in nodesVector.
							 auto [tempState, state2, pos] = *(it);
							 auto &prevNode = nodesVector[pos];
							 DDArc newArc{nextId, id, prevNode.id, decision};
							 prevNode.incomingArcs.push_back(nextId);
							 node.outgoingArcs.push_back(nextId);
							 arcs.insert(make_pair(nextId, newArc));
						 }
					 }
				 }

				 // populate nextLayer. LATER. move nodes from vector.
				 for (auto& node : nodesVector) {
					 nextLayer.push_back(node.id);
					 nodes.insert(make_pair(node.id, node));
				 }
			 }
			}
			 // build complete tree with state reduction.
			 // if (stateChangesNext) { // next DD layer will undergo state change, so create only one node.
			 //
			 // 	// if current layer is trouble maker, then handle edge case.
			 // 	const auto& layer = nodes[currentLayer.front()].globalLayer;
			 // 	bool res = stateChangesNext; //networkPtr->troubleMaker[layer] == 1;
			 //
				//  auto nextNodeId = number.getNext();
				//  DDNode newNode{nextNodeId};
				//  //newNode.nodeLayer = index;
			 //
				//  for (const auto id : currentLayer) {
				// 	 assert(nodes.count(id));
				// 	 auto &node = nodes[id];
				//  	 newNode.nodeLayer = node.nodeLayer+1;
				//  	 newNode.globalLayer = node.globalLayer+1;
				// 	 for (auto decision: node.states) {
				// 	 	if (res && node.states.size() > 1 && decision == -1) {
				// 	 		// do not create arc.
				// 	 		// cout << "Trouble maker arc is not created." << endl;
				// 	 		continue;
				// 	 	}
				// 		 auto nextId = number.getNext();
				// 		 DDArc newArc{nextId, id, nextNodeId, decision};
				// 		 node.outgoingArcs.push_back(nextId);
				// 		 newNode.incomingArcs.push_back(nextId);
				// 		 arcs.insert(make_pair(nextId, newArc));
				// 	 }
				//  } // state of this node is updated in the build() in next iteration.
				//  nextLayer.push_back(nextNodeId);
				//  nodes.insert(make_pair(nextNodeId, newNode));
			 // }
			 // else {
				//  vector<DDNode> nodesVector;
				//  unordered_set<tuple<set<int>, int, int>, tuple_hash, tuple_equal> allStates;
				//  int j = 0;
			 //
				//  for (const auto id: currentLayer) {
				// 	 assert(nodes.count(id));
				// 	 auto &node = nodes[id];
				// 	 auto statesCopy = node.states;
			 //
				// 	 for (auto decision: node.states) {
				// 		 auto newStates(statesCopy);
				// 		 if (decision != -1) newStates.erase(decision);
				// 		 auto nextId = number.getNext();
			 //
				// 		 // if newStates already in allStates, update exising node in nodesVector.
				// 		 auto [it, isInserted] = allStates.insert({newStates, 0, j});
				// 		 if (isInserted) { // create new Node
				// 			 DDNode newNode{nextId};
				// 			 DDArc newArc{nextId, id, nextId, decision};
				// 			 newNode.nodeLayer = node.nodeLayer+1;
				// 		 	 newNode.globalLayer = node.globalLayer+1;
				// 			 newNode.states = newStates;
				// 			 newNode.incomingArcs.push_back(nextId);
				// 			 node.outgoingArcs.push_back(nextId);
				// 			 arcs.insert(make_pair(nextId, newArc));
				// 			 nodesVector.push_back(newNode);
				// 			 j++;
				// 		 } else { // state already exists, update existing node in nodesVector.
				// 			 auto [tempState, state2, pos] = *(it);
				// 			 auto &prevNode = nodesVector[pos];
				// 			 DDArc newArc{nextId, id, prevNode.id, decision};
				// 			 prevNode.incomingArcs.push_back(nextId);
				// 			 node.outgoingArcs.push_back(nextId);
				// 			 arcs.insert(make_pair(nextId, newArc));
				// 		 }
				// 	 }
				//  }
			 //
				//  // populate nextLayer. LATER. move nodes from vector.
				//  for (auto& node : nodesVector) {
				// 	 nextLayer.push_back(node.id);
				// 	 nodes.insert(make_pair(node.id, node));
				//  }
			 // }
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
