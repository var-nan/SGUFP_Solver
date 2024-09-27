//
// Created by nandgate on 9/26/2024.
//

#include <vector>
#include <iostream>

#include "Network.h"
#include "gurobi_c++.h"

using namespace std;

void actualProblem(Network network) {

    // populate netowrk here
    //Network network{"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};

    int n = network.n;

    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    model.set(GRB_IntParam_InfUnbdInfo, 1);
	cout << "Environment created." << endl;
    GRBVar** X = new GRBVar * [n];
    //GRBVar** R = new GRBVar *[n];
   // vector<vector<int>> R(n); // TODO: need to populate this rewards vector.
    // for (int i = 0; i < network.networkNodes.size(); i++) {
    //     // iterate through its outgoing arcs
    //     const auto q = network.networkNodes[i];
    //     for (int k = 0; k < q.outgoingArcs.size(); k++) {
    //         int j = q.outgoingArcs[k];
    //         auto reward = network.networkArcs[j].rewards[0];
    // //        R[i][j] = reward;
    //     }
    // }
    GRBVar*** y = new GRBVar ** [n];

    for (int i = 0; i < n; i++) {
        X[i] = model.addVars(n, GRB_CONTINUOUS);
        //R[i] = model.addVars(n, GRB_CONTINUOUS);

        GRBVar** y_temp = new GRBVar *[n];

        // initialize y
        for (int j = 0; j < n; j++) {
           y_temp[j] = model.addVars(n, GRB_BINARY);
        }
        y[i] = y_temp;
    }
	cout << "Variables defined." << endl;

    model.update();

    // add objective function below.
    GRBLinExpr objFun = 0;

    for (int i = 0; i < n; i++) {
        for (int arcId : network.networkNodes[i].outgoingArcs) {
			int j = network.networkArcs[arcId].headId;
        	int reward = network.networkArcs[arcId].rewards[0];
            objFun += reward * X[i][j];
        }
    }
	model.setObjective(objFun, GRB_MAXIMIZE);
    model.update();
	cout << "Objective function defined" << endl;
    //GRBVar** beta = new GRBVar * [n];


    for (int q = 0; q < n; q++) {
		if (q==0 || q==42) {
			continue;
		}
    	GRBLinExpr c1 = 0;
		for (int incomingArc : network.networkNodes[q].incomingArcs) {
			int i = network.networkArcs[incomingArc].tailId;
			c1 += X[i][q];
		}
    	for (int outArc : network.networkNodes[q].outgoingArcs) {
			int j = network.networkArcs[outArc].headId;
			c1 -= X[q][j];
		}
    	model.addConstr(c1 == 0, "2b");
    }

	cout << "First constraint defined" << endl;

 int s = 0;
	for (int q = 0; q < n; q++){
		for (uint inArcID : network.networkNodes[q].incomingArcs)
		{
			//GRBLinExpr LHS = 0;
			uint i = network.networkArcs[inArcID].tailId;
			uint l_iq = network.networkArcs[inArcID].lowerCapacities[s];
			uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
			model.addConstr(l_iq  <= X[i][q] , "7b");
			// GRBLinExpr LHS = 0;
			model.addConstr(X[i][q] <= u_iq , "7b");
		}
	}

for (uint q : network.Vbar) {
	for (uint inArcID : network.networkNodes[q].incomingArcs) {
		for (uint outArcID : network.networkNodes[q].outgoingArcs)
		{
			uint i = network.networkArcs[inArcID].tailId;
			uint j = network.networkArcs[outArcID].headId;
			uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
			GRBLinExpr LHS = 0;
			LHS+=X[i][q]-X[q][j] + u_iq * y[i][q][j];
			model.addConstr(LHS <= u_iq , "7b");
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
			LHS+=X[q][j]-X[i][q] + u_qj * y[i][q][j];

			model.addConstr(LHS <= u_qj , "7b");
		}
	}
}
//

for (uint q : network.Vbar)
{
	for (uint inArcID : network.networkNodes[q].incomingArcs)
	{
		GRBLinExpr LHS = 0;
		uint i = network.networkArcs[inArcID].tailId;
		LHS += X[i][q];
		uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
		for (uint outArcID : network.networkNodes[q].outgoingArcs)
		{
			uint j = network.networkArcs[inArcID].headId;
			LHS -= y[i][q][j] * u_iq;
		}
		model.addConstr(LHS <= 0, "7b");
	}
}
//
for (uint q : network.Vbar)
{
	for (uint outArcID : network.networkNodes[q].outgoingArcs)
	{
		GRBLinExpr LHS = 0;
		uint j = network.networkArcs[outArcID].headId;
		LHS += X[q][j];
		uint u_q_j = network.networkArcs[outArcID].upperCapacities[s];
		for (uint inArcID : network.networkNodes[q].incomingArcs)
		{
			uint i = network.networkArcs[inArcID].tailId;
			LHS -= y[i][q][j] * u_q_j ;
		}
		model.addConstr(LHS <= 0, "7b");
	}
}

	model.update();
	model.optimize();

}
