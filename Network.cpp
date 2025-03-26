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

		uint nNodes,m, scenarios;
		file >> nNodes >> m >> scenarios;

		std::vector<NetworkNode> netNodes(nNodes);
		std::vector<NetworkArc> netArcs(m);

		for (uint i = 0; i < nNodes; i++)
			netNodes[i].nodeId = i;

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
			netNodes[tailId].outNodeIds.push_back(headId);
			netNodes[headId].inNodeIds.push_back(tailId);
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
		a4 = {};

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

		shuffleVBarNodes();

		// INFO change the order of nodes in the Vbar. Make sure that arcs of a particular node are together.
		int i = 0;
		for (const auto& id: Vbar){
			// another way of doing state update map.
			auto& node = networkNodes[id];
			set<int> states (node.outgoingArcs.begin(), node.outgoingArcs.end());
			if (node.outgoingArcs.size() < node.incomingArcs.size())
				states.insert(-1);
			stateUpdateMap.insert({i, states});
			bool first = true;
			// sort incoming arcs based in decreasing order of reward
			// std::sort(node.incomingArcs.begin(), node.incomingArcs.end(),
			// 	[this](const auto a, const auto b) constexpr {
			// 		return this->networkArcs[a].rewards[0] > this->networkArcs[b].rewards[0];
			//
			// });
			for (const auto& inId: node.incomingArcs){
				if (first) {this->hasStateChanged.push_back(true); first = false;}
				else this->hasStateChanged.push_back(false);
				processingOrder.emplace_back(i++, inId);
				layerRewards.push_back(networkArcs[inId].rewards[0]);
			}

		}
		hasStateChanged.push_back(false); // extra element for buildNextLayer in DD class.
	}
	else {
		// error in opening file. exit program.
		cerr << "File could not be opened!" << endl;
		cerr << "Error code: " << strerror(errno);
	}

}

#include <algorithm>
void Network::shuffleVBarNodes() noexcept {

	auto copyOfVbar = Vbar;
	vector<uint> demandPoints;
	vector<uint> neworder = {};

	for (auto node: networkNodes) {
		if (node.outgoingArcs.size() == 1 && networkArcs[node.outgoingArcs[0]].headId == n-1 ) {
			demandPoints.push_back(node.nodeId);
		}
	}
	for (auto id: demandPoints) {
		auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
		if (it != copyOfVbar.end()) {
			neworder.push_back(*it);
		}
	}
	for (auto id: neworder) {
		auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
		if (it != copyOfVbar.end()) {
			copyOfVbar.erase(it);
		}
	}

	vector<uint> childsVec = demandPoints;
	vector<uint> parentsVec = {};

	while (copyOfVbar.size() > 0) {
		parentsVec={};
		for (auto id: childsVec) {
			for (auto id: networkNodes[id].incomingArcs) {
				parentsVec.push_back(networkArcs[id].tailId);
			}
		}

		for (auto id: parentsVec) {
			auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
			if (it != copyOfVbar.end()) {
				auto it1 = find(neworder.begin(), neworder.end(), id);
				if (it1 == neworder.end()) {
					neworder.push_back(*it);
				}
			}
		}

		for (auto id: neworder) {
			auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
			if (it != copyOfVbar.end()) {
				copyOfVbar.erase(it);
			}
		}
		childsVec = parentsVec;
	}
	Vbar = neworder;
}