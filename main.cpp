#include <iostream>
#include "Network.h"
#include "grb.h"
#include "DD.h"
#include <chrono>
#include "OriginalProblem.h"

#include "newCGO.h"

using namespace std;

void printNetwork(const Network& network){

	cout << "Number of nodes in the network: " << network.n << "\n"
		<< "Number of edges in the network: " << network.edges << "\n"
		<< "Number of nodes in v_bar: " << network.Vbar.size() << "\n";

	const auto& networkArc = network.networkArcs[2];

	cout << networkArc.arcId << " " << networkArc.tailId << " " << networkArc.headId << " "<< networkArc.rewards.size() << endl;

	/*for (const auto& e: network.Vbar)
		cout << e << " ";
	cout << endl; */

	//testWorking();

}


size_t getNumNodesDD(const DD& dd) {
	size_t count = 0;

	for (const auto& layer : dd.tree) {
		count += layer.size();
	}
	return count;
}

void print_tree(const DD& dd) {
	for (const auto& layer: dd.tree) {
		for (const auto id: layer) cout << id << " ";
		cout << endl;
	}
}

int old_main() {
	return 0;
//
// 	/* read input and build the core datastructures */
// 	/* assuming the input is text file with first line containing two numbers n (number of nodes) and number of edges */
// 	cout << "Starting program" << endl;
// 	string fileName = "C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt";
// 	cout << fileName << endl;
// 	Network network (fileName);
//
// 	// SolveOriginalProblem(network);
// 	//
// 	// return 0;
//
// 	// auto vb = {37,22,7,30};
// 	// for (auto id: vb) {
// 	// 	for( auto incoming: network.networkNodes[22].incomingArcs) cout << incoming << " " ; cout << endl;
// 	// 	for( auto outgoing: network.networkNodes[22].outgoingArcs) cout << outgoing << " " ; cout << endl;
// 	// }
// 	//
// 	// 	cout << "Number of nodes in isNodeInVBar array: " << network.isNodeInVbar.size() << endl;
// 	//
// 	// 	int vbarcount = 0;
// 	// 	for (const auto& netNode : network.networkNodes) {
// 	// 		if (netNode.isVbar) vbarcount++;
// 	// 	}
// 	// 	assert(vbarcount == network.Vbar.size());
// 	//
// 	// 	cout << "check 1" << endl;
//
//
// 		//  // print network.
// 		// cout << "Number of nodes: " << network.n << endl;
// 		// cout << "Number of edges: " << network.edges << endl;
// 		// cout << "Number of scenarios: " << network.nScenarios << endl;
// 		// cout << endl;
// 		// for (auto& node: network.networkNodes) {
// 		// 	cout << "Node: " << node.nodeId << endl;
// 		// 	cout << "Incoming arcs: ";
// 		// 	for (auto incomingArc : node.incomingArcs) cout << incomingArc << " , "; cout << endl;
// 		// 	cout << "Outgoing arcss: "; for (auto outArc : node.outgoingArcs) cout << outArc << " , "; cout << endl;
// 		// }
// 		//
// 		// cout << "Printing arcs: "<< endl;
// 		// for (auto arc: network.networkArcs) {
// 		// 	cout << arc.arcId << " : " << arc.tailId << "->" << arc.headId << endl;
// 		// }
//
// 		//cout << "Nodes in Vbar: size: " << network.isNodeInVbar.size() << " " << network.Vbar.size() << endl;
// 		//for (auto b: network.isNodeInVbar) cout << b << " "; cout << endl;
//
// 		// return 0;
//
// 		DD restrictedDD{RESTRICTED};
// 		DDNode root{0};
// 		restrictedDD.build(root);
// 		cout << "Building tree completed oh yeah" << endl;
//
// 		vector<vector<int>> solutions = {
// 			{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
// 			{-1, -1, 71, -1, -1, -1, -1, -1, -1 ,-1},
// 			{71, -1, -1, 58, -1, 67, 10, 11, 20, 30},
// 			{71, -1, -1, 58, -1, 67, 10, 27, 30, 20},
// 			{71, -1, -1 ,58, -1, 67, 10, 20, 30, 11},
// 			{71, -1, -1, 58, 67, -1, 10, 11, 30, 20},
// 			 {-1, -1, 71, 58, 67, -1, 10 ,11, 20, 30},
// 			{-1, -1,71,58,67,-1,11,10,27,20},
// 			{-1,-1,71,58,67,-1,10,11,30,27},
// 			{-1,-1,71,58,67,-1,10,30,11,20},
// 			{-1,71,-1,58,67,-1,10,27,11,30},
// 			{-1,-1,71,58,67,-1,10,11,27,30},
// 			{-1,-1,71,58,67,-1,10,30,11,27},
// 			{-1,-1,71,58,67,-1,10,30,11,22},
// 			{-1,-1,71,58,67,-1,10,11,27,24},
// 			{-1,-1,71,58,67,-1,10,11,22,30},
// 			{-1,-1,71,58,67,-1,10,11,30,24},
// 			{-1,-1,71,58,67,-1,10,11,30,88}
//
// 		};
//
// 		// for (auto solution : solutions) {
// 		// 	cout << "Solution: " ; for (auto s: solution) cout << s << " "; cout << endl;
// 		// 	GRBEnv env;
// 		// 	GuroSolver solver{env, static_cast<int>(network.n)};
// 		// 	auto y_bar = w2y(solution, network);
// 		// 	auto cut = solver.solveSubProblemInstance(network,y_bar, 0);
// 		// 	//auto cut = generateNewCut(network, solution);
// 		// 	if (cut.cutType == OPTIMALITY) cout << "Cut type: Optimality cut" << endl;
// 		// 	else cout << "Cut type: Feasibility cut" << endl;
// 		// 	// cout << "Cut type:" << cut.cutType << endl;
// 		// 	cout << "Cut value: " << cut.RHS << endl;
// 		// 	for (auto [k,v] : cut.cutCoeff) {
// 		// 		auto[ i,q,j] = k;
// 		// 		// cout << "i: " << i << ", q: " << q << " j: " << j << "  ; val = " << v << endl;
// 		// 	}
// 		// 	//cout << "Cut coefficients" << endl;
// 		//
// 		//
// 		// 	if (cut.cutType == OPTIMALITY) {
// 		// 		restrictedDD.applyOptimalityCutRestrictedLatest(network , cut);
// 		// 	}
// 		// 	else{
// 		// 		restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut);
// 		// 		restrictedDD.displayStats();
// 		// 	}
// 		//
// 		// }
// 		//return 0;
//
// 		// auto solution = restrictedDD.solution();
// 		// for (int i = 0; i < solution.size(); ++i) solution[i] = -1;
// 		// cout << "Solution: " ;
// 		// for (auto s : solution) cout << s << " ";
// 		// cout << endl;
// 		//GRBEnv env = GRBEnv();
// 		//  GuroSolver solver{env, static_cast<int>(network.n)};
// 		//
// 		//  // auto solution  = restrictedDD.solution();
// 		//  for (int i = 0; i < solution.size(); ++i) solution[i] = -1;
// 		//  cout << "Solution: ";
// 		//  for (auto id: solution) cout << id << " , ";
// 		//  cout << endl;
// 		//  auto y_bar = w2y(solution, network);
// 		// auto cut = solver.solveSubProblemInstance(network, y_bar, 0);
// 		//
// 		// //auto [cut, heuristic] = solveSubProblemInstanceAnotherCode(network, solution);
// 		//
// 		// cout << "Cut type:" << cut.cutType << endl;
// 		// cout << "Cut value: " << cut.RHS << endl;
// 		// cout << "Cut coefficients" << endl;
// 		//
// 		// double oldSolution = std::numeric_limits<double>::max();
// 		//
// 		// double newSolution = 0;
//
// 		vi  oldWSol;
// 		vi newWSol;
//
// 		int i = 0;
// 		while (true) {
//
// 			GRBEnv tempEnv;
// 			tempEnv.set(GRB_IntParam_OutputFlag, 1);
//
// 			GuroSolver solver{tempEnv, static_cast<int>(network.n)};
//
// 			newWSol = restrictedDD.solution();
// 			if (i == 0) {
// 				for (int p = 0; p < newWSol.size(); p++) newWSol[p] = -1;
// 				//newWSol = {-1,71,-1,58,32,67,10,27,11,20};
// 			}
// 			cout << "Solution: "; for (auto x: newWSol) cout << x << " "; cout << endl;
// 			if (oldWSol != newWSol) {
// 				oldWSol = newWSol;
// 			}
// 			else {
// 				cout << "Previous Solution returned" << endl;
// 				break;
// 			}
//
// 			auto y_bar = w2y(newWSol, network);
// 			//auto cut = generateNewCut(network, newWSol);
// 			auto cut = solver.solveSubProblemInstance(y_bar, 0);
//
// 			if (cut.cutType == OPTIMALITY) cout << "Cut type: Optimality cut" << endl;
// 			else cout << "Cut type: Feasibility cut" << endl;
// 			// cout << "Cut type:" << cut.cutType << endl;
// 			cout << "Cut value: " << cut.RHS << endl;
// 			for (auto [k,v] : cut.cutCoeff) {
// 				auto[ i,q,j] = k;
// 				//cout << "i: " << i << ", q: " << q << " j: " << j << "  ; val = " << v << endl;
// 			}
// 			//cout << "Cut coefficients" << endl;
//
//
// 			if (cut.cutType == OPTIMALITY) {
// 				restrictedDD.applyOptimalityCutRestrictedLatest(cut);
// 			}
// 			else restrictedDD.applyFeasibilityCutRestrictedLatest(cut);
//
//
// 			if (restrictedDD.nodes.size() <  3 || restrictedDD.arcs.size() < 3) {
// 				cout << "*************** DELETED COMPLETE TREE **********************" << endl;
// 				break;
// 			}
//
// 			cout << "Completed " << ++i << " th iteration." << endl;
// 			cout << "After applying cut" << endl;
// 			// restrictedDD.displayStats();
//
// 		}
// 		// auto coeff = cut.cutCoeff;
//
// 		// for (auto [k,v] : coeff) {
// 		// 	auto [i,q,j] = k;
// 		// 	cout << "i: " << i << " q: " << q << " j: " << j << " val: " << v << endl;
// 		// }
// 	cout << "Program completed." << endl;
// 	return 0;
}

#include "DDSolver.h"

int main() {

	string fileName ="/home/nandgate/CLionProjects/SGUFP_Solver/42_85_20.txt";
	Network network{fileName};


	// print v bar nodes and its outoging arcs.
	// for (int id = 0; id < network.n; ++id) {
	// 	cout << id << " : incoming= " << network.networkNodes[id].incomingArcs.size()
	// 			<< " , outgoing= " << network.networkNodes[id].outgoingArcs.size() << endl;
	// }


	#ifdef DEBUG
	cout << "A4 size: " << network.A4.size() << endl;

	for (auto id : network.Vbar) {
		auto node = network.networkNodes[id];
		// cout << "Node : " << id << " inc size: " << node.incomingArcs.size() << " , outsize : " << node.outgoingArcs.size() << endl;
	}
	#endif

	SolveOriginalProblem(network);

	cout << "Solved Original problem " <<endl;
	cout << endl;
	cout << endl;

	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	cout << "Vbar order: "; for (auto id : networkPtr->Vbar) cout << id << " "; cout << endl;
	cout << "Max Width : " << MAX_WIDTH << endl;
	cout << "DEBUG enabled, solving for a subset of V Bar nodes. Single scenario." << endl;
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	using std::chrono::seconds;
	auto t1 = high_resolution_clock::now();
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	DDSolver solver{networkPtr};
	solver.initialize();
	int n_initial_cuts = 60;
	auto cuts = solver.initializeCutsParallel(n_initial_cuts);
	cout << "Number of initial cuts: " << n_initial_cuts << ". Optimality: " << cuts.second.cuts.size() <<
		" , Feasibility: " << cuts.first.cuts.size() << endl;
	cout << "**********************************************************************************************************\n\n\n" << endl;

	solver.startSolveParallel(cuts);
	auto t2 = high_resolution_clock::now();
	// cout << "Node queue strategy: LIFO" << endl;
	auto ms_int = duration_cast<seconds>(t2-t1);
	duration<double> ms_double = t2-t1;
	std::cout <<"program took " << ms_int.count() << " seconds" << endl;

	cout << "Solver finished" << endl;
	return 0;
}
