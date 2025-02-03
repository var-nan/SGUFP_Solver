/**
 * Created by nandgate on 6/2/24
 */

#include "Network.h"
#include <fstream>
#include <cstring>
#include <string>

Network::Network(const std::string& p_fileName){

	// open file.
	std::ifstream file {p_fileName};

	if (file.is_open()){
		/* file opened for reading. */
		// read nNodes and m
		uint nNodes,m; // nNodes = number of nodes, m = number of arcs // LATER change the way n,m,scenarios reads.
		uint scenarios; // number of scenarios

		file >> nNodes >> m >> scenarios;
//		nNodes = 43;
//		m = 90; //
//		scenarios = 50;
		std::vector<NetworkNode> netNodes(nNodes);
		std::vector<NetworkArc> netArcs(m);

		for (uint i = 0; i < nNodes; i++){
			// initialize a vector with nNodes empty elements.
			netNodes[i].nodeId = i;
		}

		for (uint i = 0; i < m; i++){
			// from second line, read each line and its corresponding scenario.
			// each line contains headId, headId,
			uint tailId, headId;
			file >> tailId >> headId;

			vi lowerCapacities(scenarios);
			vi upperCapacities (scenarios);
			vi rewards (scenarios);
			for (uint j = 0; j < scenarios; j++){
				// for each scenario read lb, ub, reward
				int lb, ub, reward;
				file >> lb >> ub >> reward;

				lowerCapacities[j] = lb;
				upperCapacities[j] = ub;
				rewards[j] = reward;
			}
			netArcs[i] = {i, tailId, headId, upperCapacities, lowerCapacities, rewards};
			// update node's incoming and outgoing attributes.
			netNodes[tailId].outNodeIds.push_back(headId); // netNodes[tailId].outDegree++; useful in computing heuristic in helper function.
			netNodes[headId].inNodeIds.push_back(tailId); // netNodes[headId].inDegree++;
			netNodes[tailId].outgoingArcs.push_back(i);
			netNodes[headId].incomingArcs.push_back(i);
		}

		// read v_bar
		std::string temp; file >> temp;

		// read v_bar nodes and build network.
		vector<uint> vBarNodes;

		uint index;
		//for (size_t i = 0; i < netNodes.size() ; i++) isNodeInVbar.push_back(false); // init isNodeInVbar
		while (file >> index){
			vBarNodes.push_back(index);
			netNodes[index].isVbar = true;
			//isNodeInVbar[index] = true;
		}

		// populate a1,a2,a3,a4 for subproblem formulation.
		vui a1 = {}, a2 = {}, a3 = {}, a4 = {};

		for (const auto& arc: netArcs){
			auto i = arc.tailId;
			auto j = arc.headId;
			if (netNodes[i].incomingArcs.empty()){
				if (netNodes[j].outgoingArcs.empty()) a4.push_back(arc.arcId);
				else a1.push_back(arc.arcId);
			}
			else {
				if (netNodes[j].outgoingArcs.empty()) a2.push_back(arc.arcId);
				else a3.push_back(arc.arcId);
			}
		}
		// a4 = {};
		// { // remove arcs that cannot be matched with outgoing arcs of Vbar node.
		// 	cout << "Printing whole network" << endl;
		// 	for (auto node: networkNodes) {
		// 		cout << node.nodeId << " : incoming= " << node.incomingArcs.size()
		// 		<< " , outgoing= " << node.outgoingArcs.size() << endl;
		// 	} cout << "end of network" << endl;
		//
		// 	vui arcsToRemove;
		//
		// 	for (auto inArcId : netNodes[nNodes-1].incomingArcs) {
		// 		const auto& incomingNode = netNodes[netArcs[inArcId].tailId];
		// 		if (!incomingNode.isVbar) continue;
		// 		// find max demand and compare with capacities of incoming arcs.
		// 		int maxDemandOut = netArcs[inArcId].upperCapacities[0];
		// 		for (auto qInArcId : incomingNode.incomingArcs) {
		// 			if (netArcs[qInArcId].upperCapacities[0] < maxDemandOut) arcsToRemove.push_back(qInArcId);
		// 		}
		// 	}
		//
		// 	for (auto node : netNodes) {
		// 		// remove arcs corresponding to nodes that doesn't have path from source to sink.
		// 		if ((node.outgoingArcs.empty()) &&
		// 			!(node.nodeId == 0 || node.nodeId == nNodes-1)){
		// 			// what to do?
		// 			if (!node.incomingArcs.empty())
		// 				arcsToRemove.insert(arcsToRemove.end(), node.incomingArcs.begin(), node.incomingArcs.end());
		// 			if (!node.outgoingArcs.empty())
		// 				arcsToRemove.insert(arcsToRemove.end(), node.outgoingArcs.begin(), node.outgoingArcs.end());
		// 		}
		// 	}
		// 	cout << "Number of arcs to be removed: " << arcsToRemove.size() << endl;
		// 	// for each removing arc, update head and tail nodes.
		// 	for (auto id: arcsToRemove) {
		// 		const auto& inArc = netArcs[id];
		// 		auto& inNode = netNodes[ inArc.tailId];
		// 		auto& outNode = netNodes[inArc.headId];
		// 		auto pos = find(inNode.outgoingArcs.begin(), inNode.outgoingArcs.end(), id);
		// 		if (pos != inNode.outgoingArcs.end()) inNode.outgoingArcs.erase(pos);
		// 		auto pos2 = find(outNode.incomingArcs.begin(), outNode.incomingArcs.end(), id);
		// 		if (pos2 != outNode.incomingArcs.end()) outNode.incomingArcs.erase(pos2);
		// 	}
		// }

		//reorder vBar nodes, and do not add -1 to it.

		this->n = nNodes;
		this->edges = m;
		this->networkNodes = std::move(netNodes);
		this->networkArcs = std::move(netArcs);
		this->nScenarios = scenarios;
		this->Vbar = std::move(vBarNodes);
		this->A1 = std::move(a1);
		this->A2 = std::move(a2);
		this->A3 = std::move(a3);
		this->A4 = std::move(a4);
		//this->isNodeInVbar = nodeInVbar;

		auto troubleMakerNodes = getTroubleNodes();

		shuffleVBarNodes();
		// reduce the size of vBar. just for testing.
		//this->Vbar.erase(this->Vbar.begin()+2, this->Vbar.end());

		// INFO change the order of nodes in the Vbar. Make sure that arcs of a particular node are together.
		int i = 0;
		for (const auto& id: Vbar){
			// another way of doing state update map.
			auto& node = networkNodes[id];
			set<int> states (node.outgoingArcs.begin(), node.outgoingArcs.end());
			// uncomment it later.
			if (!(node.incomingArcs.size() ==1 && node.outgoingArcs.size() ==1) &&
				(node.outgoingArcs.size() < node.incomingArcs.size())) states.insert(-1); // add -1 to states.
			// states.insert(-1);
			stateUpdateMap.insert({i, states});
			bool first = true;
			// sort incoming arcs based in decreasing order of reward
			std::sort(node.incomingArcs.begin(), node.incomingArcs.end(),
				[this](const auto a, const auto b) constexpr {
					return this->networkArcs[a].rewards[0] > this->networkArcs[b].rewards[0];

			});
			for (const auto& inId: node.incomingArcs){
				if (first) {this->hasStateChanged.push_back(true); first = false;}
				else this->hasStateChanged.push_back(false);
				processingOrder.emplace_back(i++, inId);
				layerRewards.push_back(networkArcs[inId].rewards[0]);
				troubleMaker.push_back(0);
			}
			if (find(troubleMakerNodes.begin(), troubleMakerNodes.end(), id) != troubleMakerNodes.end())
				troubleMaker.back() = 1; // set last element to 1.

		}
		hasStateChanged.push_back(false); // extra element for buildNextLayer in DD class.
		// troubleMaker.insert(troubleMaker.begin(), 0);
		cout <<"Trouble Maker structure: "; for (auto s: troubleMaker) cout << to_string(s) <<" "; cout << endl;
		cout <<"Incoming Nodes: "; for (auto id: Vbar) cout << networkNodes[id].incomingArcs.size() << " "; cout << endl;

//		// build stateUpdate map.
//		int lastId = -1;
//		for (const auto [id,aId]: processingOrder){
//			// get arc and its tail id.
//			const auto arc = netArcs[aId];
//			auto tailId = arc.tailId;
//			// get the outgoing arcs of the tail node if state needs to be changed.
//			if (lastId != tailId){
//				// state should be changed.
//				auto state = netNodes[tailId].outgoingArcs;
//				stateUpdateMap.insert({id, state}); // TODO change to aId, if needed.
//				lastId = tailId;
//			}
//		}
	}
	else {
		// error in opening file. exit program.
		cerr << "File could not be opened!" << endl;
		cerr << "Error code: " << strerror(errno);
	}

}

#include <algorithm>
void Network::shuffleVBarNodes() noexcept {
	// shift vbar nodes that  have one incoming arc and one outgoing arc to the front.
	// int i = 0, j = 0;
	//
	// while (j < this->Vbar.size()) {
	// 	const auto& node = this->networkNodes[this->Vbar[j]];
	// 	if (node.incomingArcs.size() == 1 && node.outgoingArcs.size() == 1) {
	// 		// swap i and j.
	// 		auto temp = this->Vbar[i];
	// 		this->Vbar[i] = this->Vbar[j];
	// 		this->Vbar[j] = temp;
	// 		i++;
	// 	}
	// 	j++;
	// }

	// rank : 1.05* # incomingArcs + 1.0 outgoing Arcs
	// auto rank = [this](const uint i, const uint j) {
	// 	const auto& node_i = this->networkNodes[i];
	// 	const auto& node_j = this->networkNodes[j];
	//
	// 	return (1.05 * node_i.incomingArcs.size() + 1.0 * node_i.outgoingArcs.size()) >
	// 		(1.05 * node_j.incomingArcs.size() + 1.0 * node_j.outgoingArcs.size());
	// };

	std::sort(Vbar.begin(), Vbar.end(), [this](const uint a, const uint b) {
		const auto& node_i = this->networkNodes[a];
		const auto& node_j = this->networkNodes[b];

		return (1.05 * static_cast<double>(node_i.incomingArcs.size()) + 1.0 * static_cast<double>(node_i.outgoingArcs.size())) <
			(1.05 * static_cast<double>(node_j.incomingArcs.size()) + 1.0 * static_cast<double>(node_j.outgoingArcs.size()));
	});
}

void postProcess(){
	// ASAP: change order in 'processingOrder' and 'stateUpdateMap' variables.
}

// later; this function is not necessary? remove it.
NetworkArc Network::getArc(uint32_t i, uint32_t j) const {
	// return an arc between node i (outgoing) and node j (incoming).
	for (const auto& arc: networkArcs){
		if (arc.tailId == i && arc.headId == j)
			return arc;
	}
}


vector<uint> Network::getTroubleNodes() const noexcept {

	vector<uint> troubleMakers;
	auto terminalNode = networkNodes.back().nodeId; // last node is terminal node.
	for (const auto id: Vbar) {
		// if node has arc to terminal node, trouble maker.
		const auto& node = networkNodes[id];
		// assuming trouble maker nodes have one outgoing edge and its connected to terminal.
		if (node.outNodeIds.size() == 1 && node.outNodeIds[0] == terminalNode) {
			troubleMakers.push_back(id);
		}
	}
	cout << "Trouble Maker nodes: "; for (auto id: troubleMakers) {cout << id << ", ";} cout << endl;
	return troubleMakers;
}
