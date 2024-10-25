//
// Created by nandgate on 6/5/24.
//

#include "grb.h"


void GuroSolver::initializeVariables() {

	for (int i = 0; i < n; i++){
		beta[i] = model.addVars(n, GRB_CONTINUOUS);
		gamma[i] = model.addVars(n, GRB_CONTINUOUS);
		sigma[i] = model.addVars(n, GRB_CONTINUOUS);
		phi[i] = model.addVars(n, GRB_CONTINUOUS);

		alpha[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);

		auto** lambda_temp = new GRBVar*[n];
		auto** mu_temp = new GRBVar*[n];

		for (int q = 0; q < n; q++){
			beta[i][q].set(GRB_DoubleAttr_LB, 0);
			gamma[i][q].set(GRB_DoubleAttr_LB, 0);
			sigma[i][q].set(GRB_DoubleAttr_LB, 0);
			phi[i][q].set(GRB_DoubleAttr_LB, 0);

			lambda_temp[q] = model.addVars(n, GRB_CONTINUOUS);
			mu_temp[q] = model.addVars(n, GRB_CONTINUOUS);

			for (int j = 0; j < n; j++){
				lambda_temp[q][j].set(GRB_DoubleAttr_LB, 0);
				mu_temp[q][j].set(GRB_DoubleAttr_LB, 0);
			}
		}
		lambda[i] = lambda_temp;
		mu[i] = mu_temp;
	}
	model.update();
}

void GuroSolver::addConstraints(const Network &network, int scenario) {

	/// constraint number 1 & 2 ///
	for (uint arcID : network.A1) {
		GRBLinExpr LHS = 0;
		uint q = network.networkArcs[arcID].headId;
		uint i = network.networkArcs[arcID].tailId;
		int r_iq = network.networkArcs[arcID].rewards[scenario];
		LHS += alpha[q] - beta[i][q] + gamma[i][q];
		if (network.networkNodes[q].isVbar) {
			for (uint outArc : network.networkNodes[q].outgoingArcs) {
				uint j = network.networkArcs[outArc].headId;
				LHS += lambda[i][q][j] - mu[i][q][j];
			}
			LHS += sigma[i][q];
			model.addConstr(LHS >= r_iq, "7b");
		}else {
			model.addConstr(LHS >= r_iq, "7c");
		}
	}

	/// constraint number 3 & 4 ///
	for (uint arcID : network.A2) {
		GRBLinExpr LHS = 0;
		uint j = network.networkArcs[arcID].headId;
		uint q = network.networkArcs[arcID].tailId;
		int r_qj = network.networkArcs[arcID].rewards[scenario];
		LHS += -alpha[q] - beta[q][j] + gamma[q][j];
		if (network.networkNodes[q].isVbar) {
			for (uint outArc : network.networkNodes[q].incomingArcs) {
				uint i = network.networkArcs[outArc].tailId;
				LHS += -lambda[i][q][j] + mu[i][q][j];
			}
			LHS += phi[q][j];
			model.addConstr(LHS >= r_qj, "7b");
		}else {
			model.addConstr(LHS >= r_qj, "7e");
		}
	}

	//// constraint number 5 & 6 & 7 & 8 ////
	for (uint arcID : network.A3) {
		GRBLinExpr LHS = 0;
		uint j = network.networkArcs[arcID].headId;
		uint q = network.networkArcs[arcID].tailId;
		int r_qj = network.networkArcs[arcID].rewards[scenario];
		LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
		if (network.networkNodes[q].isVbar) {
			if (network.networkNodes[j].isVbar) {
				for (uint inArc : network.networkNodes[q].incomingArcs) {
					uint i = network.networkArcs[inArc].tailId;
					LHS += mu[i][q][j] - lambda[i][q][j];
				}
				for (uint outArc : network.networkNodes[j].outgoingArcs) {
					uint i = network.networkArcs[outArc].headId;
					LHS += lambda[q][j][i] - mu[q][j][i];
				}
				LHS += phi[q][j]+sigma[q][j];
				model.addConstr(LHS >= r_qj, "7f");
			}else {
				for (uint inArc : network.networkNodes[q].incomingArcs) {
					uint i = network.networkArcs[inArc].tailId;
					LHS += mu[i][q][j] - lambda[i][q][j];
				}
				LHS +=phi[q][j];
				model.addConstr(LHS >= r_qj, "7g");
			}
		}else {
			if (network.networkNodes[j].isVbar) {
				for (uint outArc : network.networkNodes[j].outgoingArcs) {
					uint i = network.networkArcs[outArc].headId;
					LHS += lambda[q][j][i] - mu[q][j][i];
				}
				LHS += sigma[q][j];
				model.addConstr(LHS >= r_qj, "7h");
			}else {
				model.addConstr(LHS >= r_qj, "7i");
			}
		}
	}

	//// constraint number 9 ////
	// for (auto arcID1 : network.A4) {
	// 	GRBLinExpr LHS = 0;
	// 	uint q = network.networkArcs[arcID1].headId;
	// 	uint i = network.networkArcs[arcID1].tailId;
	// 	uint r_iq = network.networkArcs[arcID1].rewards[scenario];
	// 	LHS += -beta[i][q] + gamma[i][q];
	// 	model.addConstr(LHS >= r_iq, "7j");
	// }
	model.update();
}

void GuroSolver::setObjectiveFunction(const Network &network, const vector<vector<vector<shi>>> &y, int scenario) {

	GRBLinExpr obj = 0;

	const auto& nodes = network.networkNodes;
	const auto& arcs = network.networkArcs;

	// first term
	for (int q = 0; q < n; q++){
		for (uint outArc : nodes[q].outgoingArcs){
			auto j = arcs[outArc].headId;
			int u_qj = arcs[outArc].upperCapacities[scenario];
			int l_qj = arcs[outArc].lowerCapacities[scenario];
			obj += u_qj*gamma[q][j] - l_qj*beta[q][j];
		}
	}

	// second term and third term
	for (uint q : network.Vbar){
		for (uint inArc : nodes[q].incomingArcs){
			uint i = arcs[inArc].tailId;
			int u_iq = arcs[inArc].upperCapacities[scenario];
			int sum = 0;

			for (uint outArc : nodes[q].outgoingArcs) {
				//uint i = arcs[inArc].tailId;
				uint j = arcs[outArc].headId;
				int u_qj = arcs[outArc].upperCapacities[scenario];
				obj += u_iq * (1-y[i][q][j]) * lambda[i][q][j] + u_qj*(1-y[i][q][j]) * mu[i][q][j]; // second term
				sum += y[i][q][j];
			}
			// third term here.
			obj += u_iq * sum * sigma[i][q];
		}
	}

	// fourth term
	for (uint q : network.Vbar){
		for (uint outArc : nodes[q].outgoingArcs){
			auto j = arcs[outArc].headId;
			int u_qj = arcs[outArc].upperCapacities[scenario];
			int sum = 0;
			for (auto inArc : nodes[q].incomingArcs){
				uint i = arcs[inArc].tailId;
				sum += y[i][q][j];
			}
			obj += u_qj * sum * phi[q][j];
		}
	}
	model.setObjective(obj, GRB_MINIMIZE);
	model.update();
}

Cut GuroSolver::solveSubProblem(const Network &network, const vector<vector<vector<shi>>> &y_bar) {
	// TODO: run all scenarios and return feasiblity/optimality cut.
	return solveSubProblemInstance(network, y_bar, 0);
}

Cut GuroSolver::solveSubProblemInstance(const Network &network, const vector<vector<vector<shi>>> &y_bar, int scenario) {

	initializeVariables();
	//objective function.
	setObjectiveFunction(network, y_bar, scenario);
	// constraints
	addConstraints(network,scenario);
	model.optimize();

	int status = model.get(GRB_IntAttr_Status);

	if (status == GRB_OPTIMAL){
		//////////////////////////
		/// optimality cut ///////
		/////////////////////////
		CutCoefficients Y_bar_coef;
		for (uint i0 = 0; i0 < network.processingOrder.size(); i0++){
			uint arcID = network.processingOrder[i0].second;
			uint i = network.networkArcs[arcID].tailId;
			uint q = network.networkArcs[arcID].headId;
			for (uint outArc : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[outArc].headId;
				Y_bar_coef[make_tuple(i, q, j)] = 0;
			}
		}
		double rhs = 0;

		/// first term ///
		for (uint q = 0; q < network.n; q++){
			for (uint arcID : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[arcID].headId;
				rhs += network.networkArcs[arcID].upperCapacities[scenario] * gamma[q][j].get(GRB_DoubleAttr_X);
				rhs -= network.networkArcs[arcID].lowerCapacities[scenario] * beta[q][j].get(GRB_DoubleAttr_X);
			}
		}

		/// second term ///
		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					rhs += network.networkArcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_X);
					rhs += network.networkArcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_X);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_X);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_X);
				}
			}
		}

		/// third term ///
		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint u_iq = network.networkArcs[inArcID].upperCapacities[scenario];
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					Y_bar_coef[make_tuple(i, q, j)] += u_iq * sigma[i][q].get(GRB_DoubleAttr_X);
				}
			}
		}

		/// fourth term ///
		for (uint q : network.Vbar){
			for (uint outArcID : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[outArcID].headId;
				uint u_qj = network.networkArcs[outArcID].upperCapacities[scenario];
				for (uint inArcID : network.networkNodes[q].incomingArcs){
					uint i = network.networkArcs[inArcID].tailId;
					Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_X);
				}
			}
		}

		return Cut{OPTIMALITY, rhs, Y_bar_coef};
	}
	else {
		//////////////////////////////
		/////   feasibility cut   /////
		//////////////////////////////

		CutCoefficients Y_bar_coef;
		for (uint i0 = 0; i0 < network.processingOrder.size(); i0++){
			uint arcID = network.processingOrder[i0].second;
			uint i = network.networkArcs[arcID].tailId;
			uint q = network.networkArcs[arcID].headId;
			for (uint outArc : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[outArc].headId;
				// cout << "(" << i << "," << q << "," << j << ")" << endl;
				Y_bar_coef[make_tuple(i, q, j)] = 0;
			}
		}
		double rhs = 0;

		/// first term ///
		for (uint q = 0; q < network.n; q++){
			for (uint arcID : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[arcID].headId;
				rhs += network.networkArcs[arcID].upperCapacities[scenario] * gamma[q][j].get(GRB_DoubleAttr_UnbdRay);
				rhs -= network.networkArcs[arcID].lowerCapacities[scenario] * beta[q][j].get(GRB_DoubleAttr_UnbdRay);
			}
		}

		/// second term ///
		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					rhs += network.networkArcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					rhs += network.networkArcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
				}
			}
		}

		/// third term ///
		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint u_iq = network.networkArcs[inArcID].upperCapacities[scenario];
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					Y_bar_coef[make_tuple(i, q, j)] += u_iq * sigma[i][q].get(GRB_DoubleAttr_UnbdRay);
				}
			}
		}

		/// fourth term ///
		for (uint q : network.Vbar){
			for (uint outArcID : network.networkNodes[q].outgoingArcs){
				uint j = network.networkArcs[outArcID].headId;
				uint u_qj = network.networkArcs[outArcID].upperCapacities[scenario];
				for (uint inArcID : network.networkNodes[q].incomingArcs){
					uint i = network.networkArcs[inArcID].tailId;
					Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_UnbdRay);
				}
			}
		}

		return Cut{FEASIBILITY, rhs, Y_bar_coef};
	}
}