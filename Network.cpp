//
// Created by nandgate on 6/2/24.
//

#include "Network.h"
#include <fstream>


Network readNetwork(const std::string& p_fileName){
	std::ifstream file {p_fileName};

	if (file.is_open()){
		/* file opened for reading. */
		// read n and m
		uint n,m; // n = number of nodes, m = number of arcs
		usint scenarios; // number of scenarios

		file >> n >> m >> scenarios;

		std::vector<NetworkNode> networkNodes {n};
		std::vector<NetworkArc> networkArcs{m};

		for (uint i = 0; i < n; i++){
			// initialize a vector with n empty elements.
		}

		for (uint i = 0; i < m; i++){
			// read each line and its corresponding scenario.
			// each line contains " tailId, headId, weight,
			uint tailId, headId;
			file >> tailId >> headId;

			vi lowerCapacities{scenarios};
			vi upperCapacities {scenarios};
			vi rewards {scenarios};
			for (uint j = 0; j < scenarios; j++){
				// for each scenario read lb, ub, reward
				int lb, ub, reward;
				file >> lb >> ub >> reward;

				lowerCapacities.push_back(lb);
				upperCapacities.push_back(ub);
				rewards.push_back(reward);
			}
			// store this arc in the network.
			networkArcs.emplace_back(i,tailId, headId, upperCapacities, lowerCapacities, rewards);
		}

		// read v_bar
		std::string temp;
		file >> temp;

		// TODO: read v_bar nodes and build network.
	}
	else {
		// error in opening file. exit program.
	}


}

Network::Network(const std::string p_fileName) {

	Network network = readNetwork(p_fileName);


}
