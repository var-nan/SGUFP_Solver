/**
 * Created by nandgate on 6/2/24
 */

#include "Network.h"
#include <fstream>
#include <cstring>
#include <string>
#include <cmath>
Network::Network(const std::string& p_fileName){

	// open file.
	std::ifstream file {p_fileName};
	// cout << p_fileName << endl;


	if (file.is_open()) {
		/* file opened for reading. */
		// read nNodes and m
		uint nNodes,m; // nNodes = number of nodes, m = number of arcs // LATER change the way n,m,scenarios reads.
		uint scenarios; // number of scenarios
		file >> nNodes >> m >> scenarios;
		cout << "scenarios: " <<scenarios << endl;
		//		nNodes = 43;
		//		m = 90; //
		//		scenarios = 50;
		// nNodes+=2;
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
			cout << index << " * ";
		}
		cout << endl;
		cout << "got here 1 "<< endl;
		// auto troubleMakerNodes = getTroubleNodes();
		// vui arcsToRemove={};
		// // for(auto item : troubleMakerNodes){}
		// for (auto inArcId : netNodes[nNodes-1].incomingArcs) {
		// 	const auto& incomingNode = netNodes[netArcs[inArcId].tailId];
		// 	cout <<"zart" << endl;
		// 	if (!incomingNode.isVbar) continue;
		// 	// find max demand and compare with capacities of incoming arcs.
		// 	int maxDemandOut = netArcs[inArcId].upperCapacities[0];
		// 	for (auto qInArcId : incomingNode.incomingArcs) {
		// 		if (netArcs[qInArcId].upperCapacities[0] < maxDemandOut) arcsToRemove.push_back(qInArcId);
		// 	}
		// }
		// cout << "got here 2 "<< endl;
		// for (auto node : netNodes) {
		// 	// remove arcs corresponding to nodes that doesn't have path from source to sink.
		// 	if ((node.outgoingArcs.empty()) &&
		// 		!(node.nodeId == 0 || node.nodeId == nNodes-1)){
		// 		// what to do?
		// 		if (!node.incomingArcs.empty())
		// 			arcsToRemove.insert(arcsToRemove.end(), node.incomingArcs.begin(), node.incomingArcs.end());
		// 		if (!node.outgoingArcs.empty())
		// 			arcsToRemove.insert(arcsToRemove.end(), node.outgoingArcs.begin(), node.outgoingArcs.end());
		// 		}
		// }
		// cout << "Number of arcs to be removed: " << arcsToRemove.size() << endl;
		// // for each removing arc, update head and tail nodes.
		// for (auto id: arcsToRemove) {
		// 	const auto& inArc = netArcs[id];
		// 	auto& inNode = netNodes[ inArc.tailId];
		// 	auto& outNode = netNodes[inArc.headId];
		// 	auto pos = find(inNode.outgoingArcs.begin(), inNode.outgoingArcs.end(), id);
		// 	if (pos != inNode.outgoingArcs.end()) inNode.outgoingArcs.erase(pos);
		// 	auto pos2 = find(outNode.incomingArcs.begin(), outNode.incomingArcs.end(), id);
		// 	if (pos2 != outNode.incomingArcs.end()) outNode.incomingArcs.erase(pos2);
		// }
		//
		// // auto troubleMakerNodes = getTroubleNodes();
		// bool mysw =1;
		// while (mysw == 1) {
		// 	mysw =0;
		// 	for (auto item : networkNodes) {
		// 		if (item.incomingArcs.size() == 0 || item.outgoingArcs.size() == 0) {
		// 			if (item.nodeId != 0 || item.nodeId != nNodes-1) {
		// 				mysw = 1;
		// 				for (auto arc : item.outgoingArcs) {
		// 					auto pos = std::find(networkNodes[item.nodeId].outgoingArcs.begin(), networkNodes[item.nodeId].outgoingArcs.end(), arc);
		// 					networkNodes[item.nodeId].outgoingArcs.erase(pos);
		// 					pos = std::find(networkNodes[networkArcs[arc].headId].incomingArcs.begin(), networkNodes[networkArcs[arc].headId].incomingArcs.end(), arc);
		// 					networkNodes[networkArcs[arc].headId].incomingArcs.erase(pos);
		// 				}
		// 				for (auto arc : item.incomingArcs) {
		// 					auto pos = std::find(networkNodes[item.nodeId].incomingArcs.begin(), networkNodes[item.nodeId].incomingArcs.end(), arc);
		// 					networkNodes[item.nodeId].incomingArcs.erase(pos);
		// 					pos = std::find(networkNodes[networkArcs[arc].headId].outgoingArcs.begin(), networkNodes[networkArcs[arc].headId].outgoingArcs.end(), arc);
		// 					networkNodes[networkArcs[arc].tailId].outgoingArcs.erase(pos);
		// 				}
		// 			}
		// 		}
		// 	}
		// }
		// cout << "zart" <<endl;
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



		// { // remove arcs that cannot be matched with outgoing arcs of Vbar node.
		// 	cout << "Printing whole network" << endl;
		// 	for (auto node: networkNodes) {
		// 		cout << node.nodeId << " : incoming= " << node.incomingArcs.size()
		// 		<< " , outgoing= " << node.outgoingArcs.size() << endl;
		// 	} cout << "end of network" << endl;
		// 	cout << endl;
		// 	cout << endl;
		// 	cout << "Printing the V bar Nodes" << endl;
		// 	for (auto id: Vbar) {
		// 		auto node = networkNodes[id];
		// 		cout << node.nodeId << " : incoming= " << node.incomingArcs.size()
		// 		<< " , outgoing= " << node.outgoingArcs.size() << endl;
		// 	} cout << "end of V bar nodes " << endl;
		//
		// 	cout << endl;
		// 	cout << endl;
		// 	vui arcsToRemove;
		// 	cout << "printing nodes that are connected to network terminal" << endl;
		// 	for (auto inArcId : netNodes[nNodes-1].incomingArcs) {
		// 		auto incomingNode = netNodes[netArcs[inArcId].tailId];
		// 		cout << incomingNode.nodeId << " -";
		// 		if (!incomingNode.isVbar) continue;
		// 		cout<< "* ";
		// 		// find max demand and compare with capacities of incoming arcs.
		// 		int maxDemandOut = netArcs[inArcId].upperCapacities[0];
		// 		for (auto qInArcId : incomingNode.incomingArcs) {
		// 			if (netArcs[qInArcId].upperCapacities[0] < maxDemandOut) arcsToRemove.push_back(qInArcId);
		// 		}
		// 	}
		// 	cout << endl;
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


		// vector<uint> vbaredited = {};
		// for (auto item : vBarNodes) {
		// 	if(!(netNodes[item].outgoingArcs.size() == 1 && netNodes[item].incomingArcs.size() == 1)) {
		// 		vbaredited.push_back(item);
		// 	}
		// }
		//reorder vBar nodes, and do not add -1 to it.

		this->n = nNodes;
		this->edges = m;
		this->networkNodes = std::move(netNodes);
		this->networkArcs = std::move(netArcs);
		this->nScenarios = scenarios;
		// this->nScenarios = 1;
		// this->Vbar = std::move(vbaredited);
		this->Vbar = std::move(vBarNodes);
		this->A1 = std::move(a1);
		this->A2 = std::move(a2);
		this->A3 = std::move(a3);
		this->A4 = std::move(a4);

		//this->isNodeInVbar = nodeInVbar;

		// reverse(Vbar.begin(), Vbar.end());

		shuffleVBarNodes();
		// cout << 1 <<endl;
		// cout << networkNodes[Vbar.back()].incomingArcs.size();
		// while( networkNodes[Vbar.back()].incomingArcs.size() == 1 && networkNodes[Vbar.back()].outgoingArcs.size() == 1) {
		// 	// cout <<1.5<<endl;
		// 	auto item = Vbar.back();
		// 	Vbar.insert(Vbar.begin(),item);
		// 	Vbar.pop_back();
		// }

		// auto copyOfVbar = Vbar;
		// vector<uint> fisrtLayer={};
		// uint lastNodeID = nNodes - 1;
		// for (auto item : copyOfVbar) {
		// 	for(auto arc : networkNodes[item].outgoingArcs) {
		// 		if (networkArcs[arc].headId == lastNodeID) {
		// 			fisrtLayer.push_back(item);
		// 			break;
		// 		}
		// 	}
		// }
		// for (auto item : fisrtLayer) {
		// 	copyOfVbar.erase(std::find(copyOfVbar.begin(), copyOfVbar.end(), item));
		// }
		// cout << "that thing )))))))))))))))))))))" << endl;
		// for (auto item : fisrtLayer) {
		// 	cout << item << "  *";
		// }
		// cout << endl;
		// cout << "vbar" << endl;
		// for (auto item : copyOfVbar) {
		// 	cout << item << "  *";
		// }
		// cout << "incomings " << endl;
		// for (auto item: fisrtLayer) {
		// 	for (auto arc : networkNodes[item].incomingArcs) {
		// 		cout << networkArcs[arc].tailId << " ";
		// 	}
		// 	cout << " --- " ;
		// }
		vector<uint> demandPoints={};
		for (auto node: networkNodes) {
			if (node.outgoingArcs.size() == 1 && networkArcs[node.outgoingArcs[0]].headId == n-1 ) {
				demandPoints.push_back(node.nodeId);
			}
		}




		cout << endl;
		cout <<"Trouble Maker structure: "; for (auto s: troubleMaker) cout << to_string(s) <<" "; cout << endl;
		cout <<"Incoming Nodes: "; for (auto id: Vbar) cout << networkNodes[id].incomingArcs.size() << " "; cout << endl;
		cout <<"outgoing Nodes: "; for (auto id: Vbar) cout << networkNodes[id].outgoingArcs.size() << " "; cout << endl;

		// reduce the size of vBar. just for testing.
		//this->Vbar.erase(this->Vbar.begin()+2, this->Vbar.end());

		// INFO change the order of nodes in the Vbar. Make sure that arcs of a particular node are together.
		int i = 0;

		for (const auto& id: Vbar){
			// another way of doing state update map.
			auto& node = networkNodes[id];
			set<int> states (node.outgoingArcs.begin(), node.outgoingArcs.end());
			// states.insert(-1);
			// auto it = std::find(demandPoints.begin(), demandPoints.end(), id);
			if (!(node.outgoingArcs.size() >= node.incomingArcs.size())) states.insert(-1); // add -1 to states.

			// if (it != demandPoints.end() && networkNodes[i].outgoingArcs.size() == 1) {
			// 	auto maxDemand = -10000;
			// 	for(auto item : states) {
			// 		if (item  == -1 ) continue;
			// 		for(auto item2 : networkArcs[item].upperCapacities) {
			// 			if (maxDemand < networkArcs[item].upperCapacities[item2]) {
			// 				maxDemand = networkArcs[item].upperCapacities[item2];
			// 			}
			// 		}
			// 	}
			// 	for (auto  item : networkNodes[id].incomingArcs) {
			// 		if (maxDemand > networkArcs[item].upperCapacities[0]) {
			// 			states.erase(item);
			// 			cout << "happend" <<endl;
			// 		}
			// 	}
			// }



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
				troubleMaker.push_back(0);
			}
		}
		hasStateChanged.push_back(false); // extra element for buildNextLayer in DD class.
		// troubleMaker.insert(troubleMaker.begin(), 0);

		cout << "////////////////////" << endl;
		for (auto item: processingOrder) {
			cout << item.first << " " << item.second << " ** " << networkArcs[item.second].tailId << " --> " << networkArcs[item.second].headId << endl;
		}

		cout << "////////////////////" << endl;
		cout << "////////////////////" << endl;

		for (const auto& pair : stateUpdateMap) {
			std::cout << "Key: " << pair.first << " -> Values: ";
			for (const auto& value : pair.second) {
				std::cout << value << " ";
			}
			std::cout << std::endl;
		}
		cout << "////////////////////" << endl;cout << "////////////////////" << endl;







		int layersnum=0;
		for(auto id : Vbar) {
			layersnum += networkNodes[id].incomingArcs.size();
		}
		cout << "****** total Layers: " << layersnum << " *******"<< endl;


		totalLayers = layersnum;
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


void Network::vbarReOrder(bool Change) {
	stateUpdateMap={};
	if (Change) {
		Vbar.insert(Vbar.end(), Vbar[0]);
		Vbar.erase(Vbar.begin());
	}

	hasStateChanged={};
	processingOrder={};
	layerRewards={};
	troubleMaker={};

	int i = 0;
	for (const auto& id: Vbar){
		// another way of doing state update map.
		auto& node = networkNodes[id];
		set<int> states (node.outgoingArcs.begin(), node.outgoingArcs.end());
		states.insert(-1);
		// if (!(node.incomingArcs.size() ==1 && node.outgoingArcs.size() ==1) && (node.outgoingArcs.size() < node.incomingArcs.size())) states.insert(-1); // add -1 to states.
		// if (node.outgoingArcs.size() < node.incomingArcs.size()) states.insert(-1);
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
	}
	hasStateChanged.push_back(false); // extra element for buildNextLayer in DD class.

}
#include <algorithm>
void Network::shuffleVBarNodes() noexcept {
	auto copyOfVbar = Vbar;
	vector<uint> demandPoints;
	vector<uint> neworder = {};

	for (auto id : copyOfVbar) {
		cout << id << " * ";
	}
	cout << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
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
	// cout << "copyOfVbar vec:" << endl;
	// for (auto id : copyOfVbar) {
	// 	cout << id << " * ";
	// }
	// cout << endl;
	vector<uint> childsVec = demandPoints;
	vector<uint> parentsVec = {};

	while (copyOfVbar.size() > 0) {
		parentsVec={};
		for (auto id: childsVec) {
			for (auto id: networkNodes[id].incomingArcs) {
				parentsVec.push_back(networkArcs[id].tailId);
			}
		}
		// cout << "parent vec:" << endl;
		// for (auto id : parentsVec) {
		// 	cout << id << " * ";
		// }
		// cout << endl;
		for (auto id: parentsVec) {
			auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
			if (it != copyOfVbar.end()) {
				auto it1 = find(neworder.begin(), neworder.end(), id);
				if (it1 == neworder.end()) {
					neworder.push_back(*it);
				}
			}
		}
		// cout << "neworder vec:" << endl;
		// for (auto id : neworder) {
		// 	cout << id << " * ";
		// }
		// cout << endl;
		for (auto id: neworder) {
			auto it = find(copyOfVbar.begin(), copyOfVbar.end(), id);
			if (it != copyOfVbar.end()) {
				copyOfVbar.erase(it);
			}
		}
		// cout << "copyOfVbar vec:" << endl;
		// for (auto id : copyOfVbar) {
		// 	cout << id << " * ";
		// }
		// cout << endl;
		childsVec = parentsVec;
	}
	cout << "new order: " << endl;
	for (auto id : neworder) {
		cout << id << " * ";
	}
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+++++++++++++++++++++++++++++++++++++++++++++" << endl;

	Vbar = neworder;

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
