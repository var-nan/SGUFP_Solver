//
// Created by nandgate on 6/5/24.
//

#include "grb.h"

void testWorking(){
	GRBEnv env = GRBEnv(true);
	env.start();
	GRBModel model = GRBModel(env);

	std::cout << "working so far"  << std::endl;
}

inline void postProcess() {

}

void GurobiSolver::initializeVariables() {
	// initialize all variables for the current scenario.
	//
	// alpha = vector<GRBVar>(n);

	for (int i = 0; i < n; i++){
		alpha[i] = this->model.addVar(-GRB_INFINITY, GRB_INFINITY, NULL, GRB_CONTINUOUS);

		beta[i] = vector<GRBVar>(n);
		gamma[i] = vector<GRBVar>(n);
		sigma[i] = vector<GRBVar>(n);
		phi[i] = vector<GRBVar>(n);

		lambda[i] = vector<vector<GRBVar>>(n);
		mu[i] = vector<vector<GRBVar>>(n);

		for (int j = 0; j < n; j++){
			beta[i][j] = this->model.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
			gamma[i][j] = this->model.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
			sigma[i][j] = this->model.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
			phi[i][j] = this->model.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);

			lambda[i][j] = vector<GRBVar>(n);
			mu[i][j] = vector<GRBVar>(n);

			for (int k = 0; k < n; k++){
				lambda[i][j][k] = this->model.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
				mu[i][j][k] = this->model.addVar(0,GRB_INFINITY, NULL, GRB_CONTINUOUS);
			}
		}
	}

	/*
	alpha = this->model.addVars(n, GRB_CONTINUOUS);

	// todo: initialize beta, gamma, sigma, phi variables.


	for (int i = 0; i < n; i++){
		beta[i] = this->model.addVars(n, GRB_CONTINUOUS);
		gamma[i] = this->model.addVars(n, GRB_CONTINUOUS);
		sigma[i] = this->model.addVars(n, GRB_CONTINUOUS);
		phi[i] = this->model.addVars(n, GRB_CONTINUOUS);

		// TODO: initialize lambda and mu variables

		alpha[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);

		for (int j = 0; j < n; j++){
			beta[i][j].set(GRB_DoubleAttr_LB, 0);
			gamma[i][j].set(GRB_DoubleAttr_LB, 0);
			sigma[i][j].set(GRB_DoubleAttr_LB, 0);
			phi[i][j].set(GRB_DoubleAttr_LB, 0);

			lambda[i][j] = this->model.addVars(n, GRB_CONTINUOUS);
			mu[i][j] = this->model.addVars(n, GRB_CONTINUOUS);

			for (int k = 0; k < n; k++){
				//
				lambda[i][j][k].set(GRB_DoubleAttr_LB, 0);
				mu[i][j][k].set(GRB_DoubleAttr_LB, 0);
			}
		}
	}
	 */
}

/**
 * Solves the sub-problem instance for the given scenario.
 *
 * 1. Initializes variables (reset variables if defined already)
 * 2. Builds objective function.
 * 3. Runs the solver.
 * 4.
 */
void GurobiSolver::solveSubProblem(const Network& network, uint32_t scenario) {
	/*
	 * Create a model for each
	 */

	initializeVariables();

	// TODO: add objective function
	GRBLinExpr objectiveFunc;

	for (uint32_t q = 0; q < network.n; q++){

		if (network.networkNodes[q].outDegree < 1) continue;

		// iterate over outgoing nodes
		for (uint32_t j: network.networkNodes[q].outArcIds){
			auto id = network.networkArcs[j].headId;
			// extract lowerbound and upperbound for this scenario.
			int lb = network.networkArcs[j].lowerCapacities[scenario];
			int ub = network.networkArcs[j].upperCapacities[scenario];

			objectiveFunc += (-lb * beta[q][id]) + (ub * gamma[q][id]);
		}
	}

	// remaining terms
	for (const auto& q: network.Vbar){
		/* second term in objective function */
		for (const auto& in_id: network.networkNodes[q].inArcIds){

			for (const auto& out_id: network.networkNodes[q].outArcIds){
				auto i = network.networkArcs[in_id].tailId;
				auto j = network.networkArcs[out_id].headId;
				auto ub_iq = network.networkArcs[in_id].upperCapacities[scenario];
				auto ub_qj = network.networkArcs[out_id].upperCapacities[scenario];
				objectiveFunc += ub_iq * (100 /* ASAP replace with y */) + (ub_qj * (100)*(100)); /* ASAP insert lambda, mu and note values here */
			}
		}
		/* third term in objective function */
		for (const auto& in_id: network.networkNodes[q].inArcIds){
			auto i = network.networkArcs[in_id].tailId;
			auto ub = network.networkArcs[in_id].upperCapacities[scenario];
			int temp = 0;
			for (const auto& out_id: network.networkNodes[q].outArcIds) {
				// ASAP over y variables.
			}
			objectiveFunc += temp * ub * sigma[q][i];
		}

		/* fourth term in objective function */
		for (const auto& out_id: network.networkNodes[q].outArcIds){
			auto j = network.networkArcs[out_id].headId;
			auto ub = network.networkArcs[out_id].upperCapacities[scenario];
			int temp = 0;

			for (const auto& in_id: network.networkNodes[q].inArcIds){
				// ASAP get y_variables
			}
			objectiveFunc += temp * ub * phi[q][j];
		}
	}

	// INFO objective function completed.
	this->model.setObjective(objectiveFunc, GRB_MINIMIZE);

	// INFO add constraints
	// constraint 7c
	for (const auto& arc: network.networkArcs){
		auto i = arc.tailId;
		auto j = arc.headId;

		if ((network.networkNodes[j].outDegree == 0) ||
			(network.networkNodes[i].inDegree != 0)) continue; // arc subset A_1.

		GRBLinExpr LHS = alpha[j] - beta[i][j] + gamma[i][j];
		double reward = arc.rewards[scenario];
		this->model.addConstr(LHS >= reward,"7(c)");
	}

	// INFO constraint 7b
	for (const auto& arc: network.networkArcs){
		auto i = arc.tailId;
		auto q = arc.headId;

		if ((!network.networkNodes[q].isVbar) ||(network.networkNodes[q].outDegree == 0) ||
		    (network.networkNodes[i].inDegree != 0)) continue; // arc subset A_1.

		GRBLinExpr LHS = alpha[q] - beta[i][q] + gamma[i][q] + sigma[i][q];
		for (const auto& out_id: network.networkNodes[q].outArcIds)
			LHS += lambda[i][q][out_id] - mu[i][q][out_id];

		double reward = arc.rewards[scenario];
		this->model.addConstr(LHS >= reward,"7(b)");
	}

	// info: constraint 7(d), 7(e).
	for (const auto& arc: network.networkArcs) {
		auto q = arc.tailId;
		auto j = arc.headId;
		if ((network.networkNodes[q].inDegree == 0)
			|| (network.networkNodes[j].outDegree != 0)) continue; // arc set A_2

		GRBLinExpr LHS = gamma[q][j] -alpha[q] - beta[q][j];
		double reward = arc.rewards[scenario];
		//this->model.addConstr(LHS >= reward,"7(e)");

		if (network.networkNodes[q].isVbar){
			LHS += phi[q][j];
			for (const auto& i: network.networkNodes[q].inArcIds){
				LHS += mu[i][q][j] - lambda[i][q][j];
			}
		}
		this->model.addConstr(LHS >= reward, "7(d), 7(e)");
	}

	// info constraint 7(f)
	for (const auto& arc: network.networkArcs){
		auto q = arc.tailId;
		auto j = arc.headId;
		if ((network.networkNodes[q].inDegree == 0) || (network.networkNodes[j].outDegree == 0)) continue;

		// ASAP complete this constraint.
		GRBLinExpr LHS = alpha[j] - alpha[q] - beta[q][j] + gamma[q][j];
		double reward = arc.rewards[scenario];

		if ((network.networkNodes[q].isVbar) && network.networkNodes[j].isVbar) {

		}
		else if (network.networkNodes[q].isVbar){

		}
		else if (network.networkNodes[j].isVbar){ // 7(h)
			for (const auto& arc_i: network.networkNodes[q].outArcIds){
				auto i = network.networkArcs[arc_i].headId;
				LHS += lambda[q][j][i] - mu[q][j][i];
			}
		}
	}

	// INFO end constraints here.
	this->model.update(); // update changes to the model.
	this->model.optimize(); // SOLVE MODEL.

	int status = this->model.get(GRB_IntAttr_Status);

	if (status == GRB_OPTIMAL) {
		// optimality cut.
	} else {
		// feasibility cut.
	}
}
