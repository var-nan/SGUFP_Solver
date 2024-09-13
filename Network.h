/**
 * Created by nandgate on 6/1/24
 */

#pragma once

#include <utility>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <cstdint>
#include <string>

using namespace std;

typedef uint16_t usint;
typedef uint32_t uint;
typedef uint64_t ulint;
typedef vector<int> vi;
typedef vector<uint> vui;



class NetworkArc {
public:
	uint arcId;
	uint tailId; // Id of the tail node of the network.
	uint headId; // Id of the head node of the network.
	/* upper and lower capacities and corresponding rewards for different scenarios */
	vi upperCapacities;
	vi lowerCapacities;
	vi rewards;
	/* TODO: create constructors */
	NetworkArc(){}

	NetworkArc(uint id, uint tail, uint head, vi uCap, vi lCap, vi r)
		: arcId{id}, tailId{tail}, headId{head}, upperCapacities{std::move(uCap)}, lowerCapacities{std::move(lCap)}, rewards{std::move(r)} {}

};

class NetworkNode {
public:
	uint nodeId;
	//uint inDegree; // INFO removed inDegree, outDegree, inNodeIds, outNodeIds fields.
	//uint outDegree;
	//vector<uint> inNodeIds;
	//vector<uint> outNodeIds;
	vui incomingArcs; // store the id of incoming arcs in the network.
	vui outgoingArcs; // store the id of outgoing arcs in the network.
	bool isVbar;
	/* TODO: create constructors */
	NetworkNode(){}

	NetworkNode(uint id):nodeId{id}{}

	//NetworkNode(uint id, uint inDeg, uint outDeg, vui&& inNodes, vui&& outNodes)
	//	: nodeId{id}, inDegree{inDeg}, outDegree{outDeg}, inNodeIds{inNodes}, outNodeIds{outNodes} , isVbar{false}{}

	//NetworkNode(uint id, uint inDeg, uint outDeg, vui&& inArcs, vui&& outArcs, bool isV)
	//		: nodeId{id}, inDegree{inDeg}, outDegree{outDeg}, inNodeIds{inArcs}, outNodeIds{outArcs}, isVbar{isV}{}
};

class Network {

public:
	uint32_t n;
	uint32_t edges;
	vector<NetworkNode> networkNodes;
	vector<NetworkArc> networkArcs;

	//vector<std::vector<uint>> adjacencyList{n};
	vector<uint> Vbar; /* list of vertices that are in Vbar */
	uint nScenarios;
	// arc processing order
	vector<pair<int, int>> processingOrder; // should contain list of arcs in some order.
	// each element in the processing order corresponds to a layer in the DD.
	unordered_map<int, unordered_set<int>> stateUpdateMap; // (processingorder.id, vector of states)/
	// if key is not present in the stateUpdateMap, then no need to change the state at the respective arc.


	/*Network(uint nNodes, uint nEdges, vector<NetworkNode>&& netNodes,
			vector<NetworkArc>&& netArcs, unordered_set<uint>&& v_bar, uint scenarios)
		: n{nNodes}, edges{nEdges}, networkNodes{netNodes}, networkArcs{netArcs},
		Vbar{v_bar}, nScenarios{scenarios} {} */

	explicit Network(const std::string& p_fileName);

	inline NetworkArc getArc(uint32_t i, uint32_t j) const;


	/* TODO: implement move constructor ASAP*/

};

