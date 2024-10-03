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
	vector<ulint> currentLayer; // should be base layer.

	startTree = index; // LATER, size of solution vector of the root might be appropriate

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
		//i++;
		vector<ulint> nextLayer;
		nextLayer.reserve(MAX_WIDTH);
		bool isExact = buildNextLayer(currentLayer, nextLayer, index);
		if (isExact) exactLayer++; // at last, this number should be exact layer number.
		//reduceLayer(nextLayer); // INFO not doing reduction.
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}
	// terminal node layer. current layer points to last layer
	vector<ulint> terminalLayer;
	DDNode terminalNode {number.getNext()};

	for (const auto& id: currentLayer){
		// only add single arc for each node in the last layer.
		ulint arcId = number.getNext();
		auto& parentNode = nodes[id];
		// create arc and a corresponding node.
		DDArc arc{arcId, id, terminalNode.id, 1};
		arc.weight = 999999;
		terminalNode.incomingArcs.push_back(arcId);
		parentNode.outgoingArcs.push_back(arcId);
		arcs.insert(make_pair(arcId, arc));
	}

	nodes.insert(make_pair(terminalNode.id, terminalNode));
	terminalLayer.push_back(terminalNode.id);
	tree.push_back(terminalLayer);

	// update objective value for
}

inline void DD::updateState(const vector<ulint> &currentLayer, const unordered_set<int> &states){

	for (auto id: currentLayer){
		auto& node = nodes[id];
		node.states.clear();
		node.states.insert(states.begin(), states.end());
	}
}

bool DD::buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, int index) {
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
					if (count >= MAX_WIDTH) {isExact = false; break; }

					auto lastInserted = number.getNext();
					DDNode node{lastInserted};
					DDArc arc{lastInserted, parentNode.id, node.id, decision};
					node.states = parentStates;
					if (decision != -1) node.states.erase(decision);
					//node.solutionVector = parentNode.solutionVector; // solutions are computed during bulding cutset.
					//node.solutionVector.emplace_back(decision);
					node.incomingArcs.emplace_back(arc.id);
					parentNode.outgoingArcs.push_back(arc.id);
					// insert node and arc to map
					nodes.insert(std::make_pair(node.id, node));
					arcs.insert(std::make_pair(arc.id, arc));

					count++;
					//lastInserted++;

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
					//lastInserted++;
				}
			}
		}
		#elif PRUNE == RANDOM
		// merge all outgoing arcs of parent to the same children.

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

void DD::reduceLayer(vector<ulint> &currentLayer) {
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

	// delete node if no incoming arcs for the head node. TODO; INFO confirm with erfan later.
}

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
		DDArc newArc{lastInserted, dupNode.id, childNodeId, 0};
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
	double maxVal = INT32_MIN;
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


vi DD::computePathForExactNode(ulint nodeId){
	/*
	 * retrieve solution vector for a given node. used while building cutset ndoes.
	 * TODO need to verify again.
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
	cutsetNodes.reserve(MAX_WIDTH);
	int i = 0;
	for (const auto id: tree[exactLayer]){
		auto node = nodes[id];
		node.solutionVector = computePathForExactNode(id);
		//cutsetNodes[i] = node;
		cutsetNodes.push_back(node);
	}
	return cutsetNodes;
}


static vi helperFunction(const Network &network, const Cut &cut) {
	/*
	 * generates the lower bounds for each DD layer. used in feasibilitycut refinement.
	 */
	vector<int> lowerBounds(network.processingOrder.size());
	const auto& coeff = cut.cutCoeff;

	for (const auto[id_t, arcId]: network.processingOrder){
		const auto& arc = network.networkArcs[arcId];
		auto i = arc.tailId;
		auto q = arc.headId;

		const auto& node = network.networkNodes[q];
		int max = 0;
		for (const auto j: node.outNodeIds){
			auto c = coeff.at(make_tuple(static_cast<int>(i),static_cast<int>(q),static_cast<int>(j)));
			if (c > max) max = c;
		}
		lowerBounds[id_t] = max;
	}
	// compute suffix sum
	int pref = 0;
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

	vi lowerBounds = helperFunction(network, cut);

	const auto& coeff = cut.cutCoeff;
	double RHS = cut.RHS;

}

void DD::applyOptimalityCut(const Network &network, const Cut &cut) {
	// TODO: incorporate semiroot partial solution.

}

void DD::topDownDelete(ulint id) { // hard delete function.
	/*
	 * start from the current node and recursively delete the children nodes
	 * until next node is terminal node or node with multiple incoming arcs.
	 *
	 * remove node and its corresponding incoming arc together.
	 * remove the last outgoing arc of the last node.
	 *
	 * this node must have only one incoming arc (parent).
	 */
	{
		auto &node = nodes[id];
		// remove the incoming arc here.
		deleteArcById(node.incomingArcs[0]);

		for (auto outArc: node.outgoingArcs) {
			// for each outArc, find its head and apply topDownDelete() if head has single incoming arc.
			auto childId = arcs[outArc].head;
			auto &childNode = nodes[childId];
			if (childNode.incomingArcs.size() == 1) {
				topDownDelete(childId);
			} else {
				// this childId has multiple incoming arcs, just remove this arc.
				deleteArcById(outArc);
			}
		}
	}
	// here, node has no outgoing arcs and incoming arcs left,
	// also remove the node id from the tree layer.
	deleteNodeById(id);
}

void DD::removeNode(ulint id){
	/*
	 * NOTE:
	 */
	// if current node has multiple parents and multiple children, do this.
	const auto& node = nodes[id];
	//if (node.outgoingArcs.size() > 1){ // apply hard delete on eligible outgoing nodes INFO including all nodes.
	auto outArcs = node.outgoingArcs;
	for (auto childArcId : outArcs){
		const auto& childArc = arcs[childArcId];
		const auto& child = nodes[childArc.head];
		if (child.incomingArcs.size() > 1) {
			// only delete this arc
			deleteArcById(childArcId);
		}
		else { // apply hard delete on this child
			topDownDelete(child.id);
		}
	}
	//}
	//if (node.incomingArcs.size() > 1) { // multiple parents
	auto incomingArcs = node.incomingArcs;
	for (auto arcId: incomingArcs){
		const auto& arc = arcs[arcId];
		const auto& parentNode = nodes[arc.tail];
		if (parentNode.outgoingArcs.size() > 1){
			deleteArcById(arc.id);
		}
		else { // soft delete this node.
			bottomUpDelete(parentNode.id);
		}
	}
	//}
	// actual delete.
	deleteNodeById(id);

	// remove deleted ids from the tree.
	int n_removed = deletedNodeIds.size();
	auto& deletedNodeIds_l = this->deletedNodeIds;
	auto f = [&n_removed, &deletedNodeIds_l] (ulint x) mutable {
		if (deletedNodeIds_l.count(x)) { n_removed--; return true;}
		return false;
	};

	for (int i = tree.size()-2; i > 0; i--){
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


void DD::bottomUpDelete(ulint id){
	// INFO this node might contain multiple incoming parents, but might contain one (or zero) children.

	// for each incoming arc, remove arc and call soft delete on incoming node (iff has single outgoing arc).
	const auto& node = nodes[id];
	auto incomingArcs = node.incomingArcs;
	for (const auto& arcId : incomingArcs){
		const auto& arc = arcs[arcId];
		const auto& parentNode = nodes[arc.tail];
		deleteArcById(arcId); // delete this incoming arc.
		if (parentNode.outgoingArcs.empty()){
			bottomUpDelete(parentNode.id);
		}
	}
	// delete this node from nodes.
	deleteNodeById(id);

	// phase 1: recursively reach parents if they don't have multiple children.
	// if multiple incoming arcs? all bottomUpDelete () on each of the incoming nodes.
	// phase 2: recursively reach down the tree until reaching terminal node or children with multiple incoming arcs.
	// delete each node and arc in between.
//
//	auto& curNode = nodes[id];
//	// if node has multiple children now,
//	// case: when the current node has multiple incoming arcs and multiple outgoing arcs, apply deleteArc on the other nodes and
//
//	if (curNode.incomingArcs.size() > 1){
//		// for all the other arcs, apply the delete node functoin.
//		for (size_t i = 0; i < curNode.incomingArcs.size(); i++){
//			auto arcId = curNode.incomingArcs[i];
//			const auto& arc = arcs[arcId];
//			const auto& parentNode = nodes[arc.tail];
//			if (parentNode.outgoingArcs.size() > 1){
//				// do not reach this node, just remove this arc.
//				deleteArcById(arc.id);
//			}
//			else {
//				//
//				bottomUpDelete(parentNode.id);
//			}
//		}
//		const auto firstParentArcId = curNode.incomingArcs[0];
//		const auto& firstParentArc = arcs[firstParentArcId];
//		const auto& parentNode = nodes[firstParentArc.tail];
//		if (parentNode.outgoingArcs.size() > 1){
//
//		}
//		curNode = nodes[firstParentArc.tail]; // curNode points to parent node.
//	}
//
//	if (curNode.outgoingArcs.size() > 1){ // stop here and all hard delete
//	}
//
//	//const auto& arcId = curNode.incomingArcs[0];
//	//const auto& arc = arcs[arcId];
//	//auto& parentNode = nodes[arc.tail];
//
//	while (curNode.outgoingArcs.size() == 1){
//		// recursively reach parent's parent
//		if (curNode.incomingArcs.size() > 1){
//			for (size_t i = 1; i < curNode.incomingArcs.size(); i++){
//				auto parentIncomingArcId = curNode.incomingArcs[i];
//				const auto& parentIncomingArc = arcs[parentIncomingArcId];
//				const auto& parentsParentNode = nodes[parentIncomingArc.tail];
//				if (parentsParentNode.outgoingArcs.size() > 1){
//					// don't go up
//					deleteArcById(parentIncomingArcId);
//				}
//				else {
//					bottomUpDelete(parentsParentNode.id);
//				}
//			}
//		}
//		// process first arc.
//		const auto firstArcId = curNode.incomingArcs[0];
//		const auto firstArc = arcs[firstArcId];
//		curNode = nodes[firstArc.tail];
//	}
//
//	// call hard delete on the parent.
//


//	if (curNode.incomingArcs.size() > 1){
//		// multiple parents, call the function recursively
//		for (auto arcId: curNode.incomingArcs){
//			// get arc and call delete function on its tail Node.
//			auto& arc = arcs[arcId];
//			// get the parent node and check if it has multiple children
//			const auto& parentNode = nodes[arc.tail];
//			if (nodes[arc.tail].outgoingArcs.size() > 1) continue;
//			else bottomUpDelete(arc.tail);
//		}
//	}
//	else {
//		// current node has single parent.
//		while (curNode.incomingArcs.size() == 1){
//			auto& arc = arcs[curNode.incomingArcs[0]].tail;
//
//		}
//	}


}