/**
 * Created by nandgate on 6/2/24
 */

#include "Network.h"
#include <fstream>
#include <cstring>

Network::Network(const std::string& p_fileName){

	// open file.
	std::ifstream file {p_fileName};

	if (file.is_open()){
		/* file opened for reading. */
		// read nNodes and m
		uint nNodes,m; // nNodes = number of nodes, m = number of arcs // LATER change the way n,m,scenarios reads.
		uint scenarios; // number of scenarios

		//file >> nNodes >> m >> scenarios;
		nNodes = 43;
		m = 90; //
		scenarios = 50;
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
			netNodes[tailId].outNodeIds.push_back(headId); // netNodes[tailId].outDegree++;
			netNodes[headId].inNodeIds.push_back(tailId); // netNodes[headId].inDegree++;
			netNodes[tailId].outgoingArcs.push_back(i);
			netNodes[headId].incomingArcs.push_back(i);
		}

		// read v_bar
		std::string temp; file >> temp;

		// read v_bar nodes and build network.
		vector<uint> vBarNodes;

		uint index;
		while (file >> index){
			vBarNodes.push_back(index);
			netNodes[index].isVbar = true;
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

		// reduce the size of vBar. just for testing.
		//this->Vbar.erase(this->Vbar.begin()+2, this->Vbar.end());

		// INFO change the order of nodes in the Vbar. Make sure that arcs of a particular node are together.
		int i = 0;
		for (const auto& id: Vbar){
			// another way of doing state update map.
			const auto& node = networkNodes[id];
			unordered_set<int> states (node.outgoingArcs.begin(), node.outgoingArcs.end());
			states.insert(-1); // add -1 to states.
			stateUpdateMap.insert({i, states});
			for (const auto& inId: node.incomingArcs){
				processingOrder.emplace_back(i++, inId);
			}
		}

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