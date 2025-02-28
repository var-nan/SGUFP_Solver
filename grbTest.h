//
// Created by nandgate on 10/12/2024.
//

#pragma once

#include <iostream>
#include <vector>

#include "gurobi_c++.h"
#include "Network.h"

using namespace std;

void objChangeTest() {

	auto n = 10;
    vector<int> firstOBJ={1,2,3,4,5,6,7,8,9,10};
    vector<int> secondOBJ={-1,-2,-3,-4,-5,-6,-7,-8,-9,-10};

    vector<int> Acoeffs={1,2,1,2,1,2,1,2,1,2};

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_InfUnbdInfo, 1);

	GRBVar* X = model.addVars(n, GRB_CONTINUOUS);


	model.update();

	// add objective function below.
	GRBLinExpr objFun = 0;

	for (int q = 0; q < firstOBJ.size(); q++) {

		objFun += firstOBJ[q] * X[q];

	}
	model.setObjective(objFun, GRB_MINIMIZE);
	model.update();

	GRBLinExpr LHS = 0;
	for (int q = 0; q < n; q++) {
		LHS += Acoeffs[q] * X[q];
	}
	model.addConstr(LHS <= 500, "2b");
	model.addConstr(LHS >=  100, "2b");


	model.update();
	model.optimize();

	//// new obj function ///
	objFun = 0;
	for (int q = 0; q < firstOBJ.size(); q++) {

		objFun += secondOBJ[q] * X[q];

	}
	model.setObjective(objFun, GRB_MINIMIZE);
	model.update();


	model.setObjective(objFun, GRB_MINIMIZE);
	model.update();
	model.optimize();
}

