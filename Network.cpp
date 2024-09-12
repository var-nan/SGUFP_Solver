/**
 * Created by nandgate on 6/2/24
 */

#include "Network.h"
#include <fstream>
#include <iostream>


Network::Network(const std::string& p_fileName){

	// open file.
	std::ifstream file {p_fileName};

	if (file.is_open()){
		/* file opened for reading. */
		// read nNodes and m
		uint nNodes,m; // nNodes = number of nodes, m = number of arcs
		uint scenarios; // number of scenarios

		file >> nNodes >> m >> scenarios;

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
			netNodes[tailId].outNodeIds.push_back(headId); netNodes[tailId].outDegree++;
			netNodes[headId].inNodeIds.push_back(tailId); netNodes[headId].inDegree++;
		}

		// read v_bar
		std::string temp; file >> temp;

		// TODO: read v_bar nodes and build network.
		vector<uint> vBarNodes;

		while (file.peek() != EOF){
			uint i;
			file >> i;
			vBarNodes.push_back(i);
			netNodes[i].isVbar = true;
		}

		this->n = nNodes;
		this->edges = m;
		this->networkNodes = netNodes;
		this->networkArcs = netArcs;
		this->nScenarios = scenarios;
		this->Vbar = vBarNodes;

		// build Order vector.
		int i = 0;
		for (const auto& id: Vbar){
			for (const auto& inId: networkNodes[id].inNodeIds){
				processingOrder.emplace_back(pair(i++, id));
			}
		}
	}
	else {
		// error in opening file. exit program.
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