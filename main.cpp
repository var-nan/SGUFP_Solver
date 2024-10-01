#include <iostream>
#include "Network.h"
// #include "CGO.h"
#include "newCGO.h"
#include "DD.h"
#include "oneScenarioActualProblem.h"

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
#include "gurobi_c++.h"

int gurobi_main() {

	try {
		cout << "Program started" << endl;
		// set algorithmic parameters.
		// Create an environment
		GRBEnv env = GRBEnv(true);
		env.set("LogFile", "mip1.log");
		env.set(GRB_IntParam_Threads, 1); // setting at environment level,
		env.start();

		// Create an empty model
		GRBModel model = GRBModel(env);
		//model.set(GRB_IntParam_Threads, 1);

		// Create variables
		GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
		GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
		GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

		// Set objective: maximize x + y + 2 z
		model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

		// Add constraint: x + 2 y + 3 z <= 4
		model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

		// Add constraint: x + y >= 1
		model.addConstr(x + y >= 1, "c1");

		// Optimize model
		model.optimize();

		cout << x.get(GRB_StringAttr_VarName) << " "
			 << x.get(GRB_DoubleAttr_X) << endl;
		cout << y.get(GRB_StringAttr_VarName) << " "
			 << y.get(GRB_DoubleAttr_X) << endl;
		cout << z.get(GRB_StringAttr_VarName) << " "
			 << z.get(GRB_DoubleAttr_X) << endl;

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

	} catch(GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch(...) {
		cout << "Exception during optimization" << endl;
	}

	cout << "Program completed" << endl;
	return 0;
}

int main() {

	/* read input and build the core datastructures */
	/* assuming the input is text file with first line containing two numbers n (number of nodes) and number of edges */
	cout << "Starting program" << endl;
	string fileName = "C:/Users/erfank/CLionProjects/40_50_1.txt";
	cout << fileName << endl;
	Network network (fileName);
	for (int i=0 ; i< network.processingOrder.size(); i++) {
		cout << "i: " << i << "  ---  " << "(" << network.processingOrder[i].first << " , " << network.processingOrder[i].second << ")" << endl;
		cout << network.networkArcs[network.processingOrder[i].second].tailId << endl;
	}
	actualProblem(network);
	DDNode root;
	root.states = network.stateUpdateMap[0];
	root.incomingArcs = {};
	root.outgoingArcs = {};
	root.solutionVector= {};
	// popul manually
	for (int i = 0; i < network.Vbar.size(); i++) {
		cout << "node id: "<< network.Vbar[i] << "   insize: "<< network.networkNodes[network.Vbar[i]].inDegree<< "   outsize: " <<network.networkNodes[network.Vbar[i]].outDegree << " "<<endl;
	}
	cout << "those that do not have incoming arc" <<endl;
	for (int i = 0; i < network.networkNodes.size(); i++) {
		if (network.networkNodes[i].inDegree == 0) {
			cout << network.networkNodes[i].inDegree << endl;
			cout << "////////////////////////////////////" << endl;
		}
	}

	DD exactOne;
	exactOne.build(network,root,0);
	for (auto i:exactOne.tree) {
		cout << "layer size: "<< i.size() <<endl;
	}
	vi w_sol = exactOne.solution(network);
	for (const auto i : w_sol) {
		cout << i << " / ";
	}
	cout << endl;
	cout << endl;
	Cut newCut;

	newgenerateCut(w_sol , network , newCut , 0);
	cout << "newCut.RHS " << newCut.RHS << endl;
	for (const auto i : newCut.cutCoef) {
		const auto& key = i.first;
		double value = i.second;
		std::cout << "("
		  << std::get<0>(key) << ", "
		  << std::get<1>(key) << ", "
		  << std::get<2>(key) << ") -> "
		  << value << std::endl;
	}
	refineOptCut(newCut , exactOne,network);
	w_sol = exactOne.solution(network);
	for (const auto i : w_sol) {
		cout << i << " / ";
	}
	cout << endl;
	cout << endl;
	newgenerateCut(w_sol , network , newCut , 0);
	cout << "newCut.RHS " << newCut.RHS << endl;
	for (const auto i : newCut.cutCoef) {
		const auto& key = i.first;
		double value = i.second;
		std::cout << "("
		  << std::get<0>(key) << ", "
		  << std::get<1>(key) << ", "
		  << std::get<2>(key) << ") -> "
		  << value << std::endl;
	}
	refineOptCut(newCut , exactOne,network);
	w_sol = exactOne.solution(network);
	for (const auto i : w_sol) {
		cout << i << " / ";
	}
	cout << endl;
	cout << endl;
	newgenerateCut(w_sol , network , newCut , 0);
	cout << "newCut.RHS " << newCut.RHS << endl;
	for (const auto i : newCut.cutCoef) {
		const auto& key = i.first;
		double value = i.second;
		std::cout << "("
		  << std::get<0>(key) << ", "
		  << std::get<1>(key) << ", "
		  << std::get<2>(key) << ") -> "
		  << value << std::endl;
	}
	refineOptCut(newCut , exactOne,network);
	w_sol = exactOne.solution(network);
	for (const auto i : w_sol) {
		cout << i << " / ";
	}
	cout << endl;
	cout << endl;
	newgenerateCut(w_sol , network , newCut , 0);
	cout << "newCut.RHS " << newCut.RHS << endl;
	for (const auto i : newCut.cutCoef) {
		const auto& key = i.first;
		double value = i.second;
		std::cout << "("
		  << std::get<0>(key) << ", "
		  << std::get<1>(key) << ", "
		  << std::get<2>(key) << ") -> "
		  << value << std::endl;
	}














	// Cut newCut;
	// newgenerateCut(w_sol, network , newCut , 0);
	// 	if (newCut.type == 'o') {
	// 		cout << "type: " <<  newCut.type<<endl;
	// 		cout << "RHS: " <<  newCut.RHS<<endl;
	// 		for (const auto i : newCut.cutCoef) {
	// 			const auto& key = i.first;
	// 			double value = i.second;
	// 			std::cout << "("
	// 			  << std::get<0>(key) << ", "
	// 			  << std::get<1>(key) << ", "
	// 			  << std::get<2>(key) << ") -> "
	// 			  << value << std::endl;
	// 		}
	// 	}
	// cout << "newCut " << newCut.type << endl;
	// cout << "newCut.RHS " << newCut.RHS << endl;
	// vector<Cut> chosenCuts;
	// cout << " 55555555555555555555555555555555555555555555 " << endl;
	// const auto& lastLayer = exactOne.tree[exactOne.tree.size()-2];
	// for (int q : lastLayer) {
	// 	Cut newCut;
	// 	newgenerateCut(exactOne.nodes[q].solutionVector, network , newCut , 0);
	// 	for (int i : exactOne.nodes[q].solutionVector) {
	// 		cout << i << " == ";
	// 	}
	// 	// if (newCut.type == 'o') {
	// 		chosenCuts.push_back(newCut);
	// 		cout << "type: " <<  newCut.type<<endl;
	// 		cout << "RHS: " <<  newCut.RHS<<endl;
	// 		for (const auto i : newCut.cutCoef) {
	// 			const auto& key = i.first;
	// 			double value = i.second;
	// 			std::cout << "("
	// 			  << std::get<0>(key) << ", "
	// 			  << std::get<1>(key) << ", "
	// 			  << std::get<2>(key) << ") -> "
	// 			  << value << std::endl;
	// 		}
	// 		cout << "zart ==========================================sdgfblknbsbvlkndsflbkdklfjbdlfbn" << endl;
	// 	// }
	// }
	// cout << " 55555555555555555555555555555555555555555555 " << endl;
	// cout << "Nodes in VBAR: " << endl;
	// cout << "Size of isNodeVbar " << network.isNodeInVbar.size()<< endl;
	// for (int nodeId: network.isNodeInVbar) {
	// 	cout << nodeId << " ";
	// }
	// cout << endl;
	//
	// cout << "type: " <<  newCut.type<<endl;
	// cout << "RHS: " <<  newCut.RHS<<endl;
	// for (const auto i : newCut.cutCoef) {
	// 	const auto& key = i.first;
	// 	double value = i.second;
	// 	std::cout << "("
	// 	  << std::get<0>(key) << ", "
	// 	  << std::get<1>(key) << ", "
	// 	  << std::get<2>(key) << ") -> "
	// 	  << value << std::endl;
	// }



	// for(int i = 0; i < exactOne.tree.size(); i++) {
	// 	for(int j = 0; j < exactOne.tree[i].size(); j++) {
	// 		cout << exactOne.tree[i][j] << " ";
	// 	}
	// 	cout << endl;
	// 	cout << "----------------------" << endl;
	// }
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	//
	//
	// for(int i = 0; i < exactOne.tree.size(); i++) {
	// 	for(int j = 0; j < exactOne.tree[i].size(); j++) {
	// 		for(int const& k: exactOne.nodes[ exactOne.tree[i][j]].outgoingArcs) {
	// 			 cout << "(" << exactOne.arcs[k].tail <<","<< exactOne.arcs[k].head << ") -- " ;
	// 		}
	// 	}
	// 	cout << endl;
	// 	cout << "--------------------------------------" << endl;
	// }
	//
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	// cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
	//
	//
	// for(int i = 0; i < exactOne.tree.size(); i++) {
	// 	for(int j = 0; j < exactOne.tree[i].size(); j++) {
	// 		for(int const& k: exactOne.nodes[exactOne.tree[i][j]].states) {
	// 			cout <<  k  << " - " ;
	// 		}
	// 		cout << "!!";
	// 	}
	// 	cout << endl;
	// 	cout << "--------------------------------------" << endl;
	// }



	// for (unsigned int i : network.Vbar) {
	// 	cout << i << " - ";
	// 	cout << network.networkNodes[i].incomingArcs.size() <<  "  ---  ";
	// 	cout << network.networkNodes[i].inNodeIds.size() <<  "  ---  ";
	// 	cout << network.networkNodes[i].inDegree << endl;
	// }
	// cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	// for (unsigned int i : network.Vbar) {
	// 	cout << i << " - ";
	// 	cout << network.networkNodes[i].outgoingArcs.size() <<  "  ---  ";
	// 	cout << network.networkNodes[i].outNodeIds.size() <<  "  ---  ";
	// 	cout << network.networkNodes[i].outDegree << endl;
	// }
	//
	//
	//
	// // just for now.
	// // remove 25 and add it at front
	// network.Vbar.erase(std::remove(network.Vbar.begin(), network.Vbar.end(), 25), network.Vbar.end());
	// network.Vbar.insert(network.Vbar.begin(), 25);
	//
	// //printNetwork(network);
	// RestrictedDD dd{16};
	// dd.build(network);

	cout << "Program completed." << endl;
}