//
// Created by nandgate on 6/5/24.
//

//#ifndef SGUFP_SOLVER_GRB_H
//#define SGUFP_SOLVER_GRB_H
#pragma once
#include "gurobi_c++.h"
#include "Network.h"
#include <vector>
#include "Cut.h"
#include <memory>
#include <cstring>

using namespace std;

class GuroSolver{
	const shared_ptr<Network> networkPtr;
public:
	GRBEnv environment;

	int n;
	shi *y_bar;
	// gurobi variables
	GRBVar* alpha;
	GRBVar** beta;
	GRBVar** gamma;
	GRBVar** sigma;
	GRBVar** phi;
	GRBVar*** lambda;
	GRBVar*** mu;

	GRBModel model;

	GuroSolver(const shared_ptr<Network>& networkPtr_, const GRBEnv& env_): networkPtr{networkPtr_},
							environment(env_), n{static_cast<int>(networkPtr_->n)}, model {environment} {

		// set model parameters and initialize gurobi variables.
		model.set(GRB_IntParam_OutputFlag, 0);
		model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_IntParam_InfUnbdInfo, 1);

		alpha = model.addVars(n, GRB_CONTINUOUS);
		beta = new GRBVar*[n];
		gamma = new GRBVar*[n];
		sigma = new GRBVar*[n];
		phi = new GRBVar*[n];
		lambda = new GRBVar**[n];
		mu = new GRBVar**[n];
		y_bar = new shi[n*n*n];

		initializeVariables();
		addConstraints();
	}

	void initializeVariables();

	Cut solveSubProblem(const vector<vector<vector<shi>>> &y_bar);

	std::pair<CutType, Inavap::Cut> solveSubProblem(const vector<int16_t>& path);

	Cut solveSubProblemInstance(const vector<vector<vector<shi>>> &y_bar, int scenario);

	void addConstraints();

	~GuroSolver(){
		// clear up the heap and gurobi variables.
		delete[] y_bar;
		for (int i = 0; i < n; i++){
			delete[](beta[i]);
			delete[](gamma[i]);
			delete[](sigma[i]);
			delete[](phi[i]);
			for (int j = 0; j < n; j++){
				delete[](lambda[i][j]);
				delete[](mu[i][j]);
			}
			delete[](lambda[i]);
			delete[](mu[i]);
		}
		delete[](alpha);
		delete[](beta);
		delete[](gamma);
		delete[](sigma);
		delete[](phi);
		delete[](lambda);
		delete[](mu);
	}
};

//#endif //SGUFP_SOLVER_GRB_H
