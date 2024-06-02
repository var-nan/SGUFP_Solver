/**
 * Created by nandgate on 6/1/24
 */

#pragma include once

#include <utility>
#include <vector>
#include <unordered_set>
#include <cstdint>
#include <string>

typedef uint16_t usint;
typedef uint32_t uint;
typedef uint64_t ulint;
typedef std::vector<int> vi;
typedef std::vector<uint> vui;

class NetworkArc {
public:
	uint arcId;
	uint tailId;
	uint headId;
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
	uint inDegree;
	uint outDegree;
	std::vector<uint> inArcIds;
	std::vector<uint> outArcIds;
	bool isVbar;
	/* TODO: create constructors */
	NetworkNode(){}

	NetworkNode(uint id, uint inDeg, uint outDeg, vui&& inArcs, vui&& outArcs)
		: nodeId{id}, inDegree{inDeg}, outDegree{outDeg}, inArcIds{inArcs}, outArcIds{outArcs} , isVbar{false}{}

	NetworkNode(uint id, uint inDeg, uint outDeg, vui&& inArcs, vui&& outArcs, bool isV)
			: nodeId{id}, inDegree{inDeg}, outDegree{outDeg}, inArcIds{inArcs}, outArcIds{outArcs}, isVbar{isV}{}
};

class Network {

public:
	uint32_t n;
	uint32_t edges;
	std::vector<NetworkNode> networkNodes {n};
	std::vector<NetworkArc> networkArcs;

	std::vector<std::vector<uint>> adjacencyList{n};
	std::unordered_set<uint> Vbar; /* list of vertices that are in Vbar */
	usint nScenarios;

	Network(std::string p_fileName);

};

