/**
 * Created by nandgate on 6/2/24
 */

#include "Network.h"
#include <fstream>
#include <iostream>
#include <cstring>


Network::Network(const std::string& p_fileName){

	// open file.
	std::ifstream file {p_fileName};

	if (file.is_open()){
		/* file opened for reading. */
		// read nNodes and m
		uint nNodes,m; // nNodes = number of nodes, m = number of arcs
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
			// store this arc in the network. ASAP: number of nodes are greater than the given number, causing segmentation fault.
			netArcs[i] = {i, tailId, headId, upperCapacities, lowerCapacities, rewards};
			// update network nodes. //ASAP store nodeIds of outgoing and incoming arcs instead of arcIDs
			//netNodes[tailId].outNodeIds.push_back(headId); netNodes[tailId].outDegree++;
			//netNodes[headId].inNodeIds.push_back(tailId); netNodes[headId].inDegree++;
			netNodes[tailId].outgoingArcs.push_back(i);
			netNodes[headId].incomingArcs.push_back(i);
		}

		// read v_bar
		std::string temp; file >> temp;

		// TODO: read v_bar nodes and build network.
		vector<uint> vBarNodes;

		uint index;
		while (file >> index){
			vBarNodes.push_back(index);
			netNodes[index].isVbar = true;
		}
		/// new part start ///
		for (int i = 0; i < netNodes.size(); i++) {
			for (int j = 0; j < netNodes[i].outgoingArcs.size(); j++) {
				netNodes[i].outDegree+=1;
				netNodes[i].outNodeIds.push_back(netArcs[netNodes[i].outgoingArcs[j]].headId);
			}
			for (int j = 0; j < netNodes[i].incomingArcs.size(); j++) {
				netNodes[i].inDegree+=1;
				netNodes[i].inNodeIds.push_back(netArcs[netNodes[i].incomingArcs[j]].tailId);
			}
		}
		vector<int> a1 = {};
		vector<int> a2 = {};
		vector<int> a3 = {};
		vector<int> a4 = {};
		for (int i = 0; i < netNodes.size(); i++)
		{
			if (netNodes[i].incomingArcs.size() == 0) {
				if (netNodes[i].outgoingArcs.size() == 0)
				{
					a4.push_back(i);
				}
				else {
					a1.push_back(i);
				}
			}
			else {
				if (netNodes[i].outgoingArcs.size() == 0)
				{
					a2.push_back(i);
				}
				else {
					a3.push_back(i);
				}
			}
		}
		// vector<int> isnode={};
		// for (int i = 0; i < netNodes.size(); i++) {
		// 	isnode.push_back(0);
		// }
		// for (int i : vBarNodes) {
		// 	isnode[i]=1;
		// }
		/// new part end ///

		this->n = nNodes;
		this->edges = m;
		this->networkNodes = netNodes;
		this->networkArcs = netArcs;
		this->nScenarios = scenarios;
		this->Vbar = vBarNodes;
		this->A1 =a1;
		this->A2 =a2;
		this->A3 =a3;
		this->A4 =a4;




		// reduce the size of vBar. just for testing.
		//this->Vbar.erase(this->Vbar.begin()+20, this->Vbar.end());

		// TODO: change the order of nodes in the Vbar. Make sure that arcs of a particular node are together.
		int i = 0;
		for (const auto& id: this->Vbar){
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
		vector<int> isnode={};
		for (int i = 0; i < netNodes.size(); i++) {
			isnode.push_back(0);
		}
		for (int i : this->Vbar) {
			isnode[i]=1;
		}
		this->isNodeInVbar = isnode;
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

NetworkArc Network::getArc(uint32_t i, uint32_t j) const {
	// return an arc between node i (outgoing) and node j (incoming).
	for (const auto& arc: networkArcs){
		if (arc.tailId == i && arc.headId == j)
			return arc;
	}
}