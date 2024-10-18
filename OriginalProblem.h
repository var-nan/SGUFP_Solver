//
// Created by nandgate on 10/12/2024.
//

#pragma once

#include <iostream>
#include <vector>

#include "gurobi_c++.h"
#include "Network.h"

using namespace std;

void SolveOriginalProblem(Network network) {

	auto n = network.n;
    int s = 0;

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_InfUnbdInfo, 1);

	GRBVar** X = new GRBVar * [n];
	GRBVar*** y = new GRBVar ** [n];

	for (int i = 0; i < n; i++) {
		X[i] = model.addVars(n, GRB_CONTINUOUS);

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

	for (int q = 0; q < n; q++) {
		for (int arcId : network.networkNodes[q].outgoingArcs) {
			int j = network.networkArcs[arcId].headId;
			int reward = network.networkArcs[arcId].rewards[s];
			if (reward <= 0 ) {
				reward = reward;
			}
			objFun += reward * X[q][j];
		}
	}
	model.setObjective(objFun, GRB_MAXIMIZE);
	model.update();
	//GRBVar** beta = new GRBVar * [n];


	for (int q = 0; q < n; q++) {
		if(network.networkNodes[q].outgoingArcs.size() == 0 || network.networkNodes[q].incomingArcs.size() == 0) continue;
		GRBLinExpr LHS = 0;
		for (const int inArc : network.networkNodes[q].incomingArcs) {
			const int i = network.networkArcs[inArc].tailId;
			LHS += X[i][q];
		}
		for (const int outArc : network.networkNodes[q].outgoingArcs) {
			const int j = network.networkArcs[outArc].headId;
			LHS -= X[q][j];
		}
		model.addConstr(LHS == 0, "2b");
	}

	for (int q = 0; q < n; q++){
		for (uint inArcID : network.networkNodes[q].incomingArcs){
			uint i = network.networkArcs[inArcID].tailId;
			uint l_iq = network.networkArcs[inArcID].lowerCapacities[s];
			uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
			GRBLinExpr LHS = 0;
			LHS += X[i][q];
			model.addConstr(LHS >= l_iq , "2c");
			model.addConstr(LHS <= u_iq , "2c");
		}
	}

	for (uint q : network.Vbar) {
		for (uint inArcID : network.networkNodes[q].incomingArcs) {
			for (uint outArcID : network.networkNodes[q].outgoingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint j = network.networkArcs[outArcID].headId;
				uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
				GRBLinExpr LHS = 0;
				LHS+= X[i][q] - X[q][j] + u_iq * y[i][q][j];
				model.addConstr(LHS <= u_iq , "2d");
			}
		}
	}

	for (uint q : network.Vbar) {
		for (uint inArcID : network.networkNodes[q].incomingArcs) {
			for (uint outArcID : network.networkNodes[q].outgoingArcs) {
				uint i = network.networkArcs[inArcID].tailId;
				uint j = network.networkArcs[outArcID].headId;
				uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
				GRBLinExpr LHS = 0;
				LHS += X[q][j] - X[i][q] + u_qj * y[i][q][j];
				model.addConstr(LHS <= u_qj , "2e");
			}
		}
	}

	for (uint q : network.Vbar)	{
		for (uint inArcID : network.networkNodes[q].incomingArcs){
			GRBLinExpr LHS = 0;
			uint i = network.networkArcs[inArcID].tailId;
			LHS += X[i][q];
			uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
			for (uint outArcID : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[outArcID].headId;
				LHS -=  u_iq * y[i][q][j] ;
			}
			model.addConstr(LHS <= 0, "2f");
		}
	}

	for (uint q : network.Vbar) {
		for (uint outArcID : network.networkNodes[q].outgoingArcs)
		{
			GRBLinExpr LHS = 0;
			uint j = network.networkArcs[outArcID].headId;
			LHS += X[q][j];
			uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
			for (uint inArcID : network.networkNodes[q].incomingArcs)
			{
				uint i = network.networkArcs[inArcID].tailId;
				LHS -= u_qj * y[i][q][j] ;
			}
			model.addConstr(LHS <= 0, "2g");
		}
	}

	model.update();
	model.optimize();
}

