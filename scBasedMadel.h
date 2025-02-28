

#pragma once

#include <iostream>
#include <vector>

#include "gurobi_c++.h"
#include "Network.h"

using namespace std;

double sc_SolveOriginalProblem(Network network) {

	auto n = network.n;
    double scNum = network.nScenarios * 1.0;

	cout << "network.nScenarios: " << network.nScenarios << endl;

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_InfUnbdInfo, 1);
	model.set(GRB_DoubleParam_MIPGap, 0.0);
	model.set(GRB_DoubleParam_Heuristics, 0);
	model.set(GRB_IntParam_Presolve,0);
	model.set(GRB_IntParam_Cuts,0);

	GRBVar*** X = new GRBVar ** [n];
	GRBVar*** y = new GRBVar ** [n];
	for (int i = 0; i < n ; i++) {
		GRBVar** X_temp = new GRBVar *[n];
		for (int j = 0; j < n; j++) {
			X_temp[j] = model.addVars( network.nScenarios , GRB_CONTINUOUS);
		}
		X[i] = X_temp;
	}

	for (int i = 0; i < n; i++) {
		GRBVar** y_temp = new GRBVar *[n];
		// initialize y
		for (int j = 0; j < n; j++) {
		   y_temp[j] = model.addVars(n, GRB_BINARY);
		}
		y[i] = y_temp;
	}

	model.update();

	// add objective function below.
	GRBLinExpr objFun = 0;
	for (int s = 0 ; s < network.nScenarios ; s++) {
		for (int q = 0; q < n; q++) {
			for (int arcId : network.networkNodes[q].outgoingArcs) {
				int j = network.networkArcs[arcId].headId;
				int reward = network.networkArcs[arcId].rewards[s];
				objFun += ( 1.0 * reward / scNum )* X[q][j][s] ;
			}
		}
	}
	model.setObjective(objFun, GRB_MAXIMIZE);
	model.update();



	for (uint q : network.Vbar) {
		for (uint inArcID : network.networkNodes[q].incomingArcs) {
			GRBLinExpr LHS = 0;
			int i = network.networkArcs[inArcID].tailId;
			for (uint outArcID : network.networkNodes[q].outgoingArcs) {
				int j = network.networkArcs[outArcID].headId;
				LHS += y[i][q][j];
			}
			model.addConstr(LHS <= 1, "2ap");
		}
	}
	model.update();
	for (uint q : network.Vbar) {
		for (uint outArcID : network.networkNodes[q].outgoingArcs) {
			GRBLinExpr LHS = 0;
			int j = network.networkArcs[outArcID].headId;
			for (uint inArcID : network.networkNodes[q].incomingArcs) {
				int i = network.networkArcs[inArcID].tailId;
				LHS += y[i][q][j];
			}
			model.addConstr(LHS <= 1, "2ap");
		}
	}
	model.update();
	for (int q = 0; q < n; q++) {
		if(std::find(network.Vbar.begin(), network.Vbar.end(), q) == network.Vbar.end()) {
			for (uint outArcID : network.networkNodes[q].outgoingArcs) {
				int j = network.networkArcs[outArcID].headId;
				for (uint inArcID : network.networkNodes[q].incomingArcs) {
					int i = network.networkArcs[inArcID].tailId;
					model.addConstr(y[i][q][j] <= 0, "2ap");
				}
			}
		}
	}

	for (int s = 0; s < network.nScenarios; s++) {
		for (int q = 0; q < n; q++) {
			if(network.networkNodes[q].outgoingArcs.size() == 0 || network.networkNodes[q].incomingArcs.size() == 0) continue;
			GRBLinExpr LHS = 0;
			for (const int inArc : network.networkNodes[q].incomingArcs) {
				const int i = network.networkArcs[inArc].tailId;
				LHS += X[i][q][s];
			}
			for (const int outArc : network.networkNodes[q].outgoingArcs) {
				const int j = network.networkArcs[outArc].headId;
				LHS -= X[q][j][s];
			}
			model.addConstr(LHS == 0, "2b");
		}
	}
	model.update();
	for (int s = 0; s < network.nScenarios; s++) {
		for (int q = 0; q < n; q++){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint l_iq = network.networkArcs[inArcID].lowerCapacities[s];
				uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
				model.addConstr(X[i][q][s] >= l_iq , "2c");
				model.addConstr(X[i][q][s] <= u_iq , "2c");
			}
		}

	}
	model.update();
	for (int s = 0; s < network.nScenarios; s++) {
		for (uint q : network.Vbar) {
			for (uint inArcID : network.networkNodes[q].incomingArcs) {
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint i = network.networkArcs[inArcID].tailId;
					uint j = network.networkArcs[outArcID].headId;
					uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
					GRBLinExpr LHS = 0;
					LHS+= X[i][q][s] - X[q][j][s] + u_iq * y[i][q][j];
					model.addConstr(LHS <= u_iq , "2d");
				}
			}
		}

	}
	model.update();
	for (int s = 0; s < network.nScenarios ; s++) {
		for (uint q : network.Vbar) {
			for (uint inArcID : network.networkNodes[q].incomingArcs) {
				for (uint outArcID : network.networkNodes[q].outgoingArcs) {
					uint i = network.networkArcs[inArcID].tailId;
					uint j = network.networkArcs[outArcID].headId;
					uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
					GRBLinExpr LHS = 0;
					LHS += X[q][j][s] - X[i][q][s] + u_qj * y[i][q][j];
					model.addConstr(LHS <= u_qj , "2e");
				}
			}
		}

	}
	model.update();
	for (int s = 0; s < network.nScenarios ; s++) {
		for (uint q : network.Vbar)	{
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				GRBLinExpr LHS = 0;
				uint i = network.networkArcs[inArcID].tailId;
				LHS += X[i][q][s];
				uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					LHS -=  u_iq * y[i][q][j] ;
				}
				model.addConstr(LHS <= 0, "2f");
			}
		}

	}
	model.update();
	for (int s = 0; s < network.nScenarios ; s++) {
		for (uint q : network.Vbar) {
			for (uint outArcID : network.networkNodes[q].outgoingArcs)
			{
				GRBLinExpr LHS = 0;
				uint j = network.networkArcs[outArcID].headId;
				LHS += X[q][j][s];
				uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
				for (uint inArcID : network.networkNodes[q].incomingArcs)
				{
					uint i = network.networkArcs[inArcID].tailId;
					LHS -= u_qj * y[i][q][j] ;
				}
				model.addConstr(LHS <= 0, "2g");
			}
		}
	}





	model.update();
	model.optimize();
	// for(auto q : network.Vbar) {
	// 	for (uint inArc : network.networkNodes[q].incomingArcs) {
	// 		auto i = network.networkArcs[inArc].tailId;
	// 		for (uint outArc : network.networkNodes[q].outgoingArcs) {
	// 			auto j = network.networkArcs[outArc].headId;
	// 			if (y[i][q][j].get(GRB_DoubleAttr_X) > 0 ) {
	// 				cout << "i: " << i << " q " << q <<" j: " << j  << endl;
	// 			}
	// 		}
	// 	}
	// }

	return  model.get(GRB_DoubleAttr_ObjVal) ;
}

