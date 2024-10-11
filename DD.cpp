//
// Created by nandgate on 6/3/24.
//


#include "DD.h"

/**
 * Compiles the decision diagram tree with given node as the root.
 */
void DD::build(const Network& network, DDNode& node, int index) {
	// node parameter should be initialized before calling this function. node should contain states.
	// the given node parameter will be inserted as the root of the tree.

	const auto& arcOrder = network.processingOrder;
	const auto& stateUpdateMap = network.stateUpdateMap;

	// set the node to the root.
	vector<ulint> currentLayer; // should be root layer.

	startTree = index; // LATER, size of solution vector of the root might be appropriate

	node.incomingArcs.clear();
	node.outgoingArcs.clear();
	node.id = 0; // make new node if necessary.
	// TODO clear other node attributes.
	nodes.insert(std::make_pair(node.id, node));
	currentLayer.push_back(node.id);
	// insert root layer to tree.
	tree.push_back(currentLayer);

	auto start = arcOrder.begin() + index;
	auto end = arcOrder.end();

	for (; start < end; start++){

		auto[a,b] = *start;
		if (stateUpdateMap.count(a)) { // update state of each node in the layer.
			const auto& newStates = stateUpdateMap.at(a);
			updateState(currentLayer, newStates);
		}
		//const unordered_set<int> temp = stateUpdateMap.at(a);
		vector<ulint> nextLayer;
		nextLayer.reserve(MAX_WIDTH);
		bool isExact = buildNextLayer(currentLayer, nextLayer, index++);
		if (isExact) exactLayer++; // at last, this number should be exact layer number.
		//reduceLayer(nextLayer); // INFO not doing reduction.
		tree.push_back(nextLayer);
		currentLayer = std::move(nextLayer);
	}
	// terminal node layer.
	vector<ulint> terminalLayer;
	DDNode terminalNode {number.getNext()};
	// current layer points to last layer of tree.
	for (const auto& id: currentLayer){
		// only add single arc for each node in the last layer.
		ulint arcId = number.getNext();
		auto& parentNode = nodes[id];
		// create arc add to incoming of terminal node.
		DDArc arc{arcId, id, terminalNode.id, 1};
		arc.weight = INT32_MAX;
		terminalNode.incomingArcs.push_back(arcId);
		parentNode.outgoingArcs.push_back(arcId);
		arcs.insert(make_pair(arcId, arc));
	}

	nodes.insert(make_pair(terminalNode.id, terminalNode));
	terminalLayer.push_back(terminalNode.id);
	tree.push_back(terminalLayer);
}

inline void DD::updateState(const vector<ulint> &currentLayer, const unordered_set<int> &states){

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
bool DD::buildNextLayer(vector<ulint> &currentLayer, vector<ulint> &nextLayer, int index) {
	/*
	 * builds next layer from the given current layer.
	 * adds new child nodes and outgoing arcs to their respective maps.
	 * updates current layer's nodes and their arcs in the map.
	 * This function doesn't do reduction.
	 */

	bool isExact = true;

	if (type == RESTRICTED) {

		#if RESTRICTED_STRATEGY == 1
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
					node.nodeLayer = index;
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
		 * RELAXED_STRATEGY (1). make relaxed arcs to point to last inserted node in
		 * the current layer.
		 */
		#if RELAXED_STRATEGY == 1
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
		#endif

	}
	return isExact;
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
	for (size_t i = 0; i < node.solutionVector.size(); i++) solutionVector[i] = node.solutionVector[i];
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
		double max = INT32_MIN;
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

		auto arcsToDelete = node.outgoingArcs;

		for (auto outArc: arcsToDelete) {
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
		for (auto nodeId: IdsToBeRemoved)
			removeNode(nodeId);
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
		for (auto id: nodesToRemove) removeNode(id);
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
 * @param network 
 * @param cut 
 */
void DD::applyOptimalityCutRestricted(const Network &network, const Cut &cut) {
	/*
	 * 
	 */

	nodes[0].state2 = cut.RHS; // TODO change this during branch and bound.

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		auto netArc = network.processingOrder[layer-1].second;
		auto i_NetNodeId = network.networkArcs[netArc].tailId;
		auto q_NetNodeId = network.networkArcs[netArc].headId;

		vulint nodesToRemoved;

		for (auto id: tree[layer]) {
			auto& node = nodes[id];
			auto inArcId = node.incomingArcs[0];
			auto& arc = arcs[inArcId];
			auto parentNode = nodes[arc.tail];

			if (parentNode.states.find(arc.decision) != parentNode.states.end()) {
				auto j_NetNodeId = (arc.decision != -1) ? network.networkArcs[arc.decision].headId  : 0;
				arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
				node.state2 = arc.weight + parentNode.state2;
			}
			else nodesToRemoved.push_back(id);
		}

		for (const auto id: nodesToRemoved) removeNode(id); // does optimality cut removes nodes?
	}
}

/**
 * 
 * @param network
 * @param cut 
 */
void DD::applyOptimalityCutRelaxed(const Network &network, const Cut &cut) {

	nodes[0].state2 = cut.RHS;

	for (size_t layer = 1; layer < tree.size()-1; layer++) {
		const auto netArcId = network.processingOrder[layer-1].second;
		const auto i_NetNodeId = network.networkArcs[netArcId].tailId;
		const auto q_NetNodeId = network.networkArcs[netArcId].headId;

		vulint nodesToRemoved;
		vulint nodesToAdded;

		for (const auto id: tree[layer]) {

			auto& node = nodes[id];
			node.state2 = 0; // INFO setting state to zero before applying cut on the node.

			if (node.incomingArcs.size() > 1) {
				for (size_t i = 1; i < nodes[id].incomingArcs.size(); i++) {
					auto inArcId = nodes[id].incomingArcs[i];
					auto& arc = arcs[inArcId];
					const auto& parentNode = nodes[arc.tail];

					if (parentNode.states.find(arc.decision) != parentNode.states.end()) {
						// need to duplicate this node
						auto newNode = duplicate(node);
						arc.head = newNode.id;
						newNode.states = parentNode.states;
						if (arc.decision != -1) newNode.states.erase(arc.decision);
						auto j_NetNodeId = (arc.decision != -1) ? network.networkArcs[arc.decision].headId : 0;
						arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
						newNode.state2 = arc.weight + parentNode.state2;
						nodesToAdded.push_back(newNode.id);
						// copying and updating the child arcs is handled by the duplicate() function.
					}
					else {
						deleteArcById(inArcId);
					}
				}
			}
			// process first arc.
			auto inArcId = node.incomingArcs[0];
			auto& arc = arcs[inArcId];
			const auto& parentNode = nodes[arc.tail];
			node.incomingArcs = {inArcId}; // reset the incoming arcs.

			if (parentNode.states.find(arc.decision) != parentNode.states.end()) {
				auto j_NetNodeId = (arc.decision != -1 ) ? network.networkArcs[arc.decision].headId  : 0;
				arc.weight = (arc.decision != -1) ? cut.cutCoeff.at(make_tuple(i_NetNodeId, q_NetNodeId, j_NetNodeId)) : 0;
				// why remove the decided state from the node.
				node.state2 = parentNode.state2 + arc.weight;
			}
			else nodesToRemoved.push_back(id);
		}
		for (const auto id: nodesToRemoved) removeNode(id); // remove nodes marked for deletion
		tree[layer].insert(tree[layer].end(), nodesToAdded.begin(), nodesToAdded.end()); // add node ids to current layer.
	}
}
