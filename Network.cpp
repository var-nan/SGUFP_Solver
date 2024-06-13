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
		// read n and m
		uint n,m; // n = number of nodes, m = number of arcs
		usint scenarios; // number of scenarios

		file >> n >> m >> scenarios;

		std::vector<NetworkNode> networkNodes(n);
		std::vector<NetworkArc> networkArcs(m);

		for (uint i = 0; i < n; i++){
			// initialize a vector with n empty elements.
			networkNodes[i].nodeId = i;
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
			networkArcs[i] = {i,tailId, headId, upperCapacities, lowerCapacities, rewards};
			// update network nodes.
			networkNodes[tailId].outArcIds.push_back(i); networkNodes[tailId].outDegree++;
			networkNodes[headId].inArcIds.push_back(i); networkNodes[headId].inDegree++;
		}

		// read v_bar
		std::string temp; file >> temp;

		// TODO: read v_bar nodes and build network.
		vector<uint> vBarNodes;

		while (file.peek() != EOF){
			uint i;
			file >> i;
			vBarNodes.push_back(i);
			networkNodes[i].isVbar = true;
		}

		this->n = n;
		this->edges = m;
		this->networkNodes = networkNodes;
		this->networkArcs = networkArcs;
		this->nScenarios = scenarios;
		this->Vbar = vBarNodes;
	}
	else {
		// error in opening file. exit program.
	}

}

