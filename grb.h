//
// Created by nandgate on 6/5/24.
//

//#ifndef SGUFP_SOLVER_GRB_H
//#define SGUFP_SOLVER_GRB_H
#pragma once
#include "gurobi_c++.h"
#include "Network.h"
#include <vector>

using namespace std;

void testWorking();

// each thread should have GurobiSolver instance and gurobi environment.
class GurobiSolver{
private:
	GRBEnv environment;
	uint32_t n ;// = 100; // ASAP: change this number to number of nodes in the network.
	/* declare dual variables. */
	vector<GRBVar> alpha;

	vector<vector<GRBVar>> beta;
	vector<vector<GRBVar>> gamma;
	vector<vector<vector<GRBVar>>> lambda;
	vector<vector<vector<GRBVar>>> mu;
	vector<vector<GRBVar>> sigma;
	vector<vector<GRBVar>> phi;
	GRBModel model;
	//vector<GRBVar> variables = vector<GRBVar>(n);

	inline void postProcess();

public:
	// default constructor.
	GurobiSolver(const GRBEnv& env, uint32_t nodes): n{nodes}, environment(env) , model{environment} {
		//this->model = GRBModel(this->environment);
		this->alpha = vector<GRBVar>(n);
		this->beta = vector<vector<GRBVar>>(n);
		this->gamma = vector<vector<GRBVar>>(n);
		this->lambda = vector<vector<vector<GRBVar>>>(n);
		this->mu = vector<vector<vector<GRBVar>>>(n);
		this->phi = vector<vector<GRBVar>>(n);
		this->sigma = vector<vector<GRBVar>>(n);
	}

	void resetModel();

	void solveSubProblem(const Network& network, const vector<vector<vector<bool>>>& y, uint32_t scenario);

	void initializeVariables();

	void resetVariables();

	~GurobiSolver(){
		// TODO: clear all the variables.
	}
};


//#endif //SGUFP_SOLVER_GRB_H
