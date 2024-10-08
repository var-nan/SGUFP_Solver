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
	/// start of new code ///
	currentLayer = tree.back();
	vector<int> nextLayer={};
	DDNode Termnial;
	Termnial.id = number.getNext();
	Termnial.nodeLayer = nodes[currentLayer[0]].nodeLayer+1;
	for (auto id : currentLayer) {
		lastInserted = number.getNext();
		DDNode parent = nodes[id];
		DDArc arc{lastInserted, parent.id, Termnial.id, 1};
		arc.weight = 9999999;
		nodes[id].outgoingArcs.push_back(lastInserted);
		Termnial.incomingArcs.push_back(lastInserted);
		arcs.insert(std::make_pair(arc.id, arc));
	}
	nodes.insert(std::make_pair(Termnial.id , Termnial));
	nextLayer.push_back(Termnial.id);
	tree.push_back(nextLayer);
	/// end of new code ///
	for (auto id : tree[tree.size()-2]) {
		nodes[id].objVal = 99999999999;
	}
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

					lastInserted = number.getNext();
					DDNode node{lastInserted};
					DDArc arc{lastInserted, parentNode.id, node.id, decision};
					node.states = parentStates;
					if (decision != -1) node.states.erase(decision);
					node.solutionVector = parentNode.solutionVector;
					node.solutionVector.emplace_back(decision);
					node.incomingArcs.emplace_back(arc.id);
					node.nodeLayer = parentNode.nodeLayer+1;
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
#endif
	}
	else if(type == RELAXED) {
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
					node.nodeLayer = nodes[currentLayer[0]].nodeLayer+1;
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
					node.nodeLayer = nodes[currentLayer[0]].nodeLayer+1;
					nodes.insert(std::make_pair(node.id, node));
					arcs.insert(std::make_pair(arc.id, arc));
					nextLayer.push_back(node.id);
				}
			}
			if (nextLayer.size() > MAX_WIDTH) isExact = false;
		}
	}

	return isExact;
}

void refineOptCut(Cut &newCut , DD &DDTree ,Network& network) {
	DDTree.nodes[DDTree.tree[0][0]].state2 = newCut.RHS; // LATER: incorporate semiroot partial solution
	for(int i0=1 ; i0<DDTree.tree.size();i0++) {
		for (auto n : DDTree.tree[i0]) {
			DDTree.nodes[n].state2 = -9999999;
		}
	}
	// cout << "zart " << endl;
	for (int i0=1 ; i0<DDTree.tree.size()-1;i0++) {
		auto theArc = network.processingOrder[i0-1].second; // LATER: add offset of semiroot node before i0.
		auto parentNetNode = network.networkArcs[theArc].tailId;
		auto middleNetNode = network.networkArcs[theArc].headId;
		for (auto DDnodeID : DDTree.tree[i0]) {
			int newState = -9999999;
			for (auto arcID : DDTree.nodes[DDnodeID].incomingArcs) {
				auto DDparent = DDTree.arcs[arcID].tail;
				// // we don't have a arc corresponding to -1 in networkArcs.
				if (DDTree.arcs[arcID].decision == -1) {
					DDTree.arcs[arcID].weight = 0; //DDTree.nodes[DDparent].state2
				}
				else {
					int decisionNode = network.networkArcs[DDTree.arcs[arcID].decision].headId;
					DDTree.arcs[arcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,decisionNode)];
				}
				// cout << "new state = " << newState << endl;
				// cout << " DDTree.arcs[arcID].weight + DDTree.nodes[DDparent].state2; " <<DDTree.arcs[arcID].weight + DDTree.nodes[DDparent].state2 << endl;
				if (newState <= DDTree.arcs[arcID].weight + DDTree.nodes[DDparent].state2 ) newState = DDTree.arcs[arcID].weight + DDTree.nodes[DDparent].state2;
			}
			DDTree.nodes[DDnodeID].state2 = newState;
		}
		// for (int DDArcID : DDTree.nodes[DDTree.tree.back()[0]].incomingArcs) {}
	}
	int termnialState = -9999999;
	for (int inArcs : DDTree.nodes[ DDTree.tree.back()[0]].incomingArcs) {
		auto newvalue = DDTree.nodes[DDTree.arcs[inArcs].tail].state2;
		if (newvalue < DDTree.arcs[inArcs].weight) {
			DDTree.arcs[inArcs].weight =newvalue ;
		}
		if (termnialState <= DDTree.arcs[inArcs].weight) termnialState = DDTree.arcs[inArcs].weight;
	}
	DDTree.nodes[DDTree.tree.back()[0]].state2 = termnialState;




	// for(int i0=0 ; i0<DDTree.tree.size();i0++) {
	// 	for (auto n : DDTree.tree[i0]) {
	// 		cout << DDTree.nodes[n].state2 << " - ";
	// 	}
	// 	cout << endl;
	// 	cout << "----------------------------------------" << endl;
	// }
}

void DD::refineFeasibilityCut(Cut &newCut , DD &DDTree ,Network& network) {
	DDTree.nodes[DDTree.tree[0][0]].state2 = newCut.RHS;
	for(int i0=1 ; i0<DDTree.tree.size();i0++) {
		for (auto n : DDTree.tree[i0]) {
			DDTree.nodes[n].state2 = 0;
		}
	}
	for(int i0=0 ; i0<DDTree.tree.size();i0++) {
		for (auto n : DDTree.tree[i0]) {
			cout << DDTree.nodes[n].nodeLayer << " - ";
		}
		cout << endl;
		cout << "----------------------------------------" << endl;
	}

	for( auto item :newCut.heuristic) {
		cout << item << "///" ;
	}
	cout << endl;

	for(int i0=1 ; i0 < tree.size()-1;i0++) {

		auto theArc = network.processingOrder[i0-1].second;
		auto parentNetNode = network.networkArcs[theArc].tailId;
		auto middleNetNode = network.networkArcs[theArc].headId;
		vector<int> nodesToBeRemoved = {};
		cout << " before size of the tree: " << i0 << "    " << tree[i0].size()<< endl;
		vector<int> nodesToBeAdded = {};
		for (auto n : tree[i0]) {
			vi allIncomings = nodes[n].incomingArcs;
			// cout << "nodes[n].incomingArcs.size():  " << nodes[n].incomingArcs.size() << endl;
			cout << "nodes[n].incomingArcs.size(): " << nodes[n].incomingArcs.size()<<endl;
			if (allIncomings.size() > 1) {
				for (auto inc = 1 ; inc < allIncomings.size() ; inc++) {
					// cout << "got in multi arc for " << endl;
					int incArcID = allIncomings[inc];
					int decisionNode = arcs[incArcID].decision;
					auto parentNodeID = arcs[incArcID].tail;;
					if (nodes[parentNodeID].states.find(decisionNode) == nodes[parentNodeID].states.end()) {
						cout << "found it" << endl;
						continue;
					}
					arcs[incArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
					// cout << "arcs[incArcID].weight :" << arcs[incArcID].weight << endl;

					int newState = DDTree.nodes[parentNodeID].state2 + arcs[incArcID].weight;
					if (newState  + newCut.heuristic[nodes[n].nodeLayer] >= 0) {
						// cout << endl;
						// cout << "making new node " << endl;
						// cout << "newState: " << newState << endl;
						// cout << "nodeLayer:   " << i0 << endl;
						// cout << endl;
						int lastInserted = number.getNext();
						DDNode newnode = nodes[n];
						newnode.states = nodes[parentNodeID].states;
						newnode.id = lastInserted;
						newnode.incomingArcs = {allIncomings[inc]};
						newnode.state2 = newState;
						if (decisionNode != -1) {
							newnode.states.erase(decisionNode);
						}
						newnode.outgoingArcs={};
						// NOTE: can be furthur imporoved buy not really making all outgoing arcs based on the known state
						// Note: Do not build infeasible arcs
						for(int out : DDTree.nodes[n].outgoingArcs) {
							auto outArc = arcs[out];
							DDArc newArc = outArc;
							int lastInserteda = number.getNext();
							newArc.id = lastInserteda;
							newArc.tail = newnode.id;
							nodes[newArc.head].incomingArcs.push_back(newArc.id);
							newnode.outgoingArcs.push_back(newArc.id);
							arcs.insert(make_pair(newArc.id,newArc));
						}
						nodes.insert(make_pair(newnode.id,newnode));
						nodesToBeAdded.push_back(newnode.id);
						// cout << " node got added "<< endl;
						// tree[newnode.nodeLayer].push_back(newnode.id);
					}
				}
				/// the first incoming arc
				int incArcID = allIncomings[0];
				int decisionNode = arcs[incArcID].decision;
				auto parentNodeID = arcs[incArcID].tail;
				nodes[n].incomingArcs = {incArcID};
				if (nodes[parentNodeID].states.find(decisionNode) == nodes[parentNodeID].states.end()) {
					nodesToBeRemoved.push_back(n);
					cout << "found it 2" << endl;
				}
				if (decisionNode != -1) {
					arcs[incArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
				}else {
					arcs[incArcID].weight = 0;
				}
				int newState = DDTree.nodes[parentNodeID].state2 + arcs[incArcID].weight;
				if (newState  + newCut.heuristic[nodes[n].nodeLayer] >= 0) {
					// nodes[n].incomingArcs = {incArcID};
					nodes[n].state2 = newState;
					nodes[n].states = nodes[parentNodeID].states;
					if (decisionNode != -1) {
						nodes[n].states.erase(decisionNode);
					}
				}else {
					// cout << "deletion happened!"<< endl;
					nodesToBeRemoved.push_back(nodes[n].id);
					// removeNode(nodes[n].id);
				}

			}else { // there is only one incomming arc
				auto inArcID = allIncomings[0];
				int decisionNode = arcs[inArcID].decision;
				// cout << "decisionNode" << decisionNode << endl;
				auto parentNodeID = arcs[inArcID].tail;
				if (decisionNode != -1) {
					arcs[inArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
					int newState = DDTree.nodes[parentNodeID].state2 + arcs[inArcID].weight;
					if (newState  + newCut.heuristic[nodes[n].nodeLayer] >= 0) {
						DDTree.nodes[n].state2 = newState;
					}else {
						nodesToBeRemoved.push_back(nodes[n].id);
					}
				}else{
					DDTree.nodes[n].state2 = DDTree.nodes[parentNodeID].state2;
					if (DDTree.nodes[n].state2 + newCut.heuristic[nodes[n].nodeLayer] < 0) {
						// cout << "deletion happened!"<< endl;
						nodesToBeRemoved.push_back(nodes[n].id);
					}
				}
			}

			// cout << "tree[i0]" << tree[i0].size() << endl;
		}

		// cout << "tree[i0-1]" << tree[i0-1].size() << endl;
		for(auto item : nodesToBeAdded) {
			tree[i0].push_back(item);
		}
		cout << " after size of the tree:     " << tree[i0].size()<< endl;
		// cout << endl;
		for (auto badNode : nodesToBeRemoved) {
			removeNode(badNode);
		}
		cout << " after size of the tree:     " << tree[i0].size()<< endl;
	}


	// cout << "got here " << endl;
	for(int i0=0 ; i0<tree.size();i0++) {
		for (auto n : tree[i0]) {
			cout << nodes[n].state2 << " - ";
		}
		cout << endl;
		cout << "----------------------------------------" << endl;
	}

}
void DD::refineOptimalityCut(Cut &newCut , DD &DDTree ,Network& network) {
	DDTree.nodes[DDTree.tree[0][0]].state2 = newCut.RHS;
	for(int i0=1 ; i0<DDTree.tree.size();i0++) {
		for (auto n : DDTree.tree[i0]) {
			DDTree.nodes[n].state2 = 0;
		}
	}
	for(int i0=0 ; i0<DDTree.tree.size();i0++) {
		for (auto n : DDTree.tree[i0]) {
			cout << DDTree.nodes[n].id << " - ";
		}
		cout << endl;
		cout << "----------------------------------------" << endl;
	}
	for(int i0=1 ; i0 < tree.size()-1;i0++) {
		auto theArc = network.processingOrder[i0-1].second;
		auto parentNetNode = network.networkArcs[theArc].tailId;
		auto middleNetNode = network.networkArcs[theArc].headId;
		vector<int> nodesToBeRemoved = {};
		vector<int> nodesToBeAdded = {};
		for (auto n : tree[i0]) {
			vi allIncomings = nodes[n].incomingArcs;
			if (allIncomings.size() > 1) {
				for (auto inc = 1 ; inc < allIncomings.size() ; inc++) {
					int incArcID = allIncomings[inc];
					int decisionNode = arcs[incArcID].decision;
					auto parentNodeID = arcs[incArcID].tail;;
					if (nodes[parentNodeID].states.find(decisionNode) == nodes[parentNodeID].states.end()) {
						// cout << " this isssssss happpenning "<< endl;
						deleteArcById(incArcID);
						continue;
					}
					int lastInserted = number.getNext();
					DDNode newnode = nodes[n];
					newnode.id = lastInserted;
					newnode.states = nodes[parentNodeID].states;
					if (decisionNode != -1) {
						newnode.states.erase(decisionNode);
						arcs[incArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
					}else {
						arcs[incArcID].weight = 0;
						newnode.states = nodes[parentNodeID].states;
					}
					int newState = DDTree.nodes[parentNodeID].state2 + arcs[incArcID].weight;
					newnode.incomingArcs = {allIncomings[inc]};
					newnode.state2 = newState;
					newnode.outgoingArcs={};
					// NOTE: can be furthur imporoved buy not really making all outgoing arcs based on the known state
					// Note: Do not build infeasible arcs
					for(int out : DDTree.nodes[n].outgoingArcs) {
						auto outArc = arcs[out];
						DDArc newArc = outArc;
						int lastInserteda = number.getNext();
						newArc.id = lastInserteda;
						newArc.tail = newnode.id;
						nodes[newArc.head].incomingArcs.push_back(newArc.id);
						newnode.outgoingArcs.push_back(newArc.id);
						arcs.insert(make_pair(newArc.id,newArc));
					}
					nodes.insert(make_pair(newnode.id,newnode));
					nodesToBeAdded.push_back(newnode.id);
				}
				/// the first incoming arc
				int incArcID = allIncomings[0];
				int decisionNode = arcs[incArcID].decision;
				auto parentNodeID = arcs[incArcID].tail;
				nodes[n].incomingArcs = {incArcID};
				if (nodes[parentNodeID].states.find(decisionNode) == nodes[parentNodeID].states.end()) {
					nodesToBeRemoved.push_back(n);
				}else {
					if (decisionNode != -1) {
						arcs[incArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
						nodes[n].states.erase(decisionNode);
					}else {
						arcs[incArcID].weight = 0;
						nodes[n].states = nodes[parentNodeID].states;
					}
					int newState = DDTree.nodes[parentNodeID].state2 + arcs[incArcID].weight;
					nodes[n].state2 = newState;
				}
			}else { // there is only one incomming arc
				auto inArcID = allIncomings[0];
				int decisionNode = arcs[inArcID].decision;
				auto parentNodeID = arcs[inArcID].tail;
				if (nodes[parentNodeID].states.find(decisionNode) == nodes[parentNodeID].states.end()) {
					nodesToBeRemoved.push_back(n);
					cout << "this happend "<< endl;
				}else {
					if (decisionNode != -1) {
						arcs[inArcID].weight = newCut.cutCoef[make_tuple(parentNetNode,middleNetNode,network.networkArcs[decisionNode].headId)];
						int newState = DDTree.nodes[parentNodeID].state2 + arcs[inArcID].weight;
						DDTree.nodes[n].state2 = newState;
					}else{
						arcs[inArcID].weight=0;
						DDTree.nodes[n].state2 = DDTree.nodes[parentNodeID].state2;
					}
				}
			}
		}
		for(auto item : nodesToBeAdded) {
			tree[i0].push_back(item);
		}
		for (auto badNode : nodesToBeRemoved) {
			removeNode(badNode);
			cout<< endl;
		}
	}
	int termnialState = -9999999;
	for (int inArcs : DDTree.nodes[ DDTree.tree.back()[0]].incomingArcs) {
		auto newvalue = nodes[arcs[inArcs].tail].state2;
		if (newvalue < arcs[inArcs].weight) {
			arcs[inArcs].weight =newvalue ;
		}
		if (termnialState <= arcs[inArcs].weight) termnialState = DDTree.arcs[inArcs].weight;
	}
	nodes[tree.back()[0]].state2 = termnialState;

	for(int i0=0 ; i0<tree.size();i0++) {
		for (auto n : tree[i0]) {
			cout << nodes[n].state2  << " - ";
		}
		cout << endl;
		cout << "----------------------------------------" << endl;
	}


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
	deletedNodeIds.insert(id);
}

inline DDNode DD::duplicate(const DDNode& node){
	// clone the node. for every outgoing arc, create new arc and point it to child node.
	lastInserted = number.getNext();
	DDNode dupNode(lastInserted++);
	dupNode.state2 = node.state2;
	dupNode.states = node.states;
	// ASAP how solution vector is copied
	// incoming arc is updated by the caller.

	for (const auto& outArcId: node.outgoingArcs){
		const auto& childNodeId = arcs[outArcId].head;
		lastInserted = number.getNext();
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
		// bool feasible = true; // TODO: based on the cut, build if node is feasible.
		// if (feasible){
		DDNode newNode = duplicate(node); // set incoming arc.
		// update incoming arc and make the arc to point to new node.
		newNode.incomingArcs.push_back(arcId);
		nodes.insert(std::make_pair(newNode.id, newNode));
		arc.head = newNode.id;
		// }
		// else {
		// 	// delete arc.
		// 	deleteArcById(arcId);
		// 	auto& tailOutArcs = nodes[arc.tail].outgoingArcs;
		// 	tailOutArcs.erase(std::find(tailOutArcs.begin(), tailOutArcs.end(), arcId));
		// 	// remove arc
		// 	arcs.erase(arcId);
		// }
	}
	// remove the incoming arcs of original node.
	// compute state and keep it if node is feasible.
	// if (true){ // feasible, delete all incoming arcs except the first.
	// 	node.incomingArcs.erase(node.incomingArcs.begin()+1, node.incomingArcs.end());
	// }
	// else {
	// 	// delete node.
	// 	deleteArcById(node.incomingArcs[0]);
	// 	deleteNodeById(node.id);
	// }
}

static void printTree(const vector<vector<int>>& tree){

	for (const auto& layer: tree){
		for (const int node: layer){
			cout << node << " ";
		}
		cout << endl;
	}
}

vi DD::solution(Network network) {
	int maxVal = -999999;
	int maxNodeId = 0;
	vi path(tree.size()-2);
	int terminalNodeId = tree[tree.size()-1][0];
	const auto& terminalNode = nodes[terminalNodeId];

	for (const auto arcId: terminalNode.incomingArcs){ // find maxValue path.
		const auto& arc = arcs[arcId];
		if (arc.weight > maxVal){
			maxVal = arc.weight;
			maxNodeId = arc.tail;
		}
	}
	// cout << "tree.size(): " << tree.size() << endl;
	// maxVal and maxNodeId is populated.
	for (size_t l = tree.size()-2; l > 0; l--){ // check indexes later.
		const auto& node = nodes[maxNodeId];
		for(const auto& arcId: node.incomingArcs) {
			auto& arc = arcs[arcId];
			auto& tailNode = nodes[arc.tail];
			if (tailNode.state2 + arc.weight == node.state2) {
				// this is the optimal parent.
				maxNodeId = arc.tail;
				path[l-1] = arc.decision;
				// cout << "arc.decision" << arc.decision << endl;
				break;
			}
		}
	}
	const auto& rootNode = nodes[tree[0][0]];
	cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	for(auto i:path) {
		cout << i << " - ";
	}
	cout << endl;
	cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
	// cout << "rootNode.solutionVector.size()" <<rootNode.solutionVector.size() << endl;
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


vi DD::computeExactNodePartialSolutionVector(int nodeId){
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
		node.solutionVector = computeExactNodePartialSolutionVector(id);
		//cutsetNodes[i] = node;
		cutsetNodes.push_back(node);
	}
	return cutsetNodes;
}


void DD::topDownDelete(ulint id) { // hard delete function.

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

	deleteNodeById(id);
}
void DD::removeNode(ulint id){

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
	std::for_each(deletedNodeIds.begin(), deletedNodeIds.end(), [&num](ulint x) mutable{num.setNext(static_cast<int>(x));});
	std::for_each(deletedNodeIds.begin(), deletedNodeIds.end(), [](ulint x) {cout << x << " ";});
	deletedNodeIds.clear();


}
void DD::bottomUpDelete(ulint id){
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
	// phase
}