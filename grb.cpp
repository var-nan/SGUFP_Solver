//
// Created by nandgate on 6/5/24.
//

#include "grb.h"

#include <iomanip>

#include "NodeExplorer.h"


void GuroSolver::initializeVariables() {

	// for (int i = 0; i < n; i++){
	// 	beta[i] = model.addVars(n, GRB_CONTINUOUS);
	// 	gamma[i] = model.addVars(n, GRB_CONTINUOUS);
	// 	sigma[i] = model.addVars(n, GRB_CONTINUOUS);
	// 	phi[i] = model.addVars(n, GRB_CONTINUOUS);
	//
	// 	alpha[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
	//
	// 	auto** lambda_temp = new GRBVar*[n];
	// 	auto** mu_temp = new GRBVar*[n];
	//
	// 	for (int q = 0; q < n; q++){
	// 		beta[i][q].set(GRB_DoubleAttr_LB, 0);
	// 		gamma[i][q].set(GRB_DoubleAttr_LB, 0);
	// 		sigma[i][q].set(GRB_DoubleAttr_LB, 0);
	// 		phi[i][q].set(GRB_DoubleAttr_LB, 0);
	//
	// 		lambda_temp[q] = model.addVars(n, GRB_CONTINUOUS);
	// 		mu_temp[q] = model.addVars(n, GRB_CONTINUOUS);
	//
	// 		for (int j = 0; j < n; j++){
	// 			lambda_temp[q][j].set(GRB_DoubleAttr_LB, 0);
	// 			mu_temp[q][j].set(GRB_DoubleAttr_LB, 0);
	// 		}
	// 	}
	// 	lambda[i] = lambda_temp;
	// 	mu[i] = mu_temp;
	// }
	// model.update();
}

// void GuroSolver::addConstraints(const Network &network, int scenario) {

	// /// constraint number 1 & 2 ///
	// for (uint arcID : network.A1) {
	// 	GRBLinExpr LHS = 0;
	// 	uint q = network.networkArcs[arcID].headId;
	// 	uint i = network.networkArcs[arcID].tailId;
	// 	int r_iq = network.networkArcs[arcID].rewards[scenario];
	// 	LHS += alpha[q] - beta[i][q] + gamma[i][q];
	// 	if (network.networkNodes[q].isVbar) {
	// 		for (uint outArc : network.networkNodes[q].outgoingArcs) {
	// 			uint j = network.networkArcs[outArc].headId;
	// 			LHS += lambda[i][q][j] - mu[i][q][j];
	// 		}
	// 		LHS += sigma[i][q];
	// 		model.addConstr(LHS >= r_iq, "7b");
	// 	}else {
	// 		model.addConstr(LHS >= r_iq, "7c");
	// 	}
	// }
	//
	// /// constraint number 3 & 4 ///
	// for (uint arcID : network.A2) {
	// 	GRBLinExpr LHS = 0;
	// 	uint j = network.networkArcs[arcID].headId;
	// 	uint q = network.networkArcs[arcID].tailId;
	// 	int r_qj = network.networkArcs[arcID].rewards[scenario];
	// 	LHS += -alpha[q] - beta[q][j] + gamma[q][j];
	// 	if (network.isNodeInVbar[q]) {
	// 		for (uint outArc : network.networkNodes[q].incomingArcs) {
	// 			uint i = network.networkArcs[outArc].tailId;
	// 			LHS += -lambda[i][q][j] + mu[i][q][j];
	// 		}
	// 		LHS += phi[q][j];
	// 		model.addConstr(LHS >= r_qj, "7b");
	// 	}else {
	// 		model.addConstr(LHS >= r_qj, "7e");
	// 	}
	// }
	//
	// //// constraint number 5 & 6 & 7 & 8 ////
	// for (uint arcID : network.A3) {
	// 	GRBLinExpr LHS = 0;
	// 	uint j = network.networkArcs[arcID].headId;
	// 	uint q = network.networkArcs[arcID].tailId;
	// 	int r_qj = network.networkArcs[arcID].rewards[scenario];
	// 	LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
	// 	if (network.isNodeInVbar[q]) {
	// 		if (network.isNodeInVbar[j]) {
	// 			for (uint inArc : network.networkNodes[q].incomingArcs) {
	// 				uint i = network.networkArcs[inArc].tailId;
	// 				LHS += mu[i][q][j] - lambda[i][q][j];
	// 			}
	// 			for (uint outArc : network.networkNodes[j].outgoingArcs) {
	// 				uint i = network.networkArcs[outArc].headId;
	// 				LHS += lambda[q][j][i] - mu[q][j][i];
	// 			}
	// 			LHS += phi[q][j]+sigma[q][j];
	// 			model.addConstr(LHS >= r_qj, "7f");
	// 		}else {
	// 			for (uint inArc : network.networkNodes[q].incomingArcs) {
	// 				uint i = network.networkArcs[inArc].tailId;
	// 				LHS += mu[i][q][j] - lambda[i][q][j];
	// 			}
	// 			LHS +=phi[q][j];
	// 			model.addConstr(LHS >= r_qj, "7g");
	// 		}
	// 	}else {
	// 		if (network.isNodeInVbar[j]) {
	// 			for (uint outArc : network.networkNodes[j].outgoingArcs) {
	// 				uint i = network.networkArcs[outArc].headId;
	// 				LHS += lambda[q][j][i] - mu[q][j][i];
	// 			}
	// 			LHS += sigma[q][j];
	// 			model.addConstr(LHS >= r_qj, "7h");
	// 		}else {
	// 			model.addConstr(LHS >= r_qj, "7i");
	// 		}
	// 	}
	// }
	//
	// //// constraint number 9 ////
	// // for (auto arcID1 : network.A4) {
	// // 	GRBLinExpr LHS = 0;
	// // 	uint q = network.networkArcs[arcID1].headId;
	// // 	uint i = network.networkArcs[arcID1].tailId;
	// // 	uint r_iq = network.networkArcs[arcID1].rewards[scenario];
	// // 	LHS += -beta[i][q] + gamma[i][q];
	// // 	model.addConstr(LHS >= r_iq, "7j");
	// // }
	// model.update();
// }

// void GuroSolver::setObjectiveFunction(const Network &network, const vector<vector<vector<shi>>> &y, int scenario) {

	// GRBLinExpr obj = 0;
	//
	// const auto& nodes = network.networkNodes;
	// const auto& arcs = network.networkArcs;
	//
	// // first term
	// for (int q = 0; q < n; q++){
	// 	for (uint outArc : nodes[q].outgoingArcs){
	// 		auto j = arcs[outArc].headId;
	// 		int u_qj = arcs[outArc].upperCapacities[scenario];
	// 		int l_qj = arcs[outArc].lowerCapacities[scenario];
	// 		obj += u_qj*gamma[q][j] - l_qj*beta[q][j];
	// 	}
	// }
	//
	// // second term and third term
	// for (uint q : network.Vbar){
	// 	for (uint inArc : nodes[q].incomingArcs){
	// 		uint i = arcs[inArc].tailId;
	// 		int u_iq = arcs[inArc].upperCapacities[scenario];
	// 		int sum = 0;
	//
	// 		for (uint outArc : nodes[q].outgoingArcs) {
	// 			//uint i = arcs[inArc].tailId;
	// 			uint j = arcs[outArc].headId;
	// 			int u_qj = arcs[outArc].upperCapacities[scenario];
	// 			obj += u_iq * (1-y[i][q][j]) * lambda[i][q][j] + u_qj*(1-y[i][q][j]) * mu[i][q][j]; // second term
	// 			sum += y[i][q][j];
	// 		}
	// 		// third term here.
	// 		obj += u_iq * sum * sigma[i][q];
	// 	}
	// }
	//
	// // fourth term
	// for (uint q : network.Vbar){
	// 	for (uint outArc : nodes[q].outgoingArcs){
	// 		auto j = arcs[outArc].headId;
	// 		int u_qj = arcs[outArc].upperCapacities[scenario];
	// 		int sum = 0;
	// 		for (auto inArc : nodes[q].incomingArcs){
	// 			uint i = arcs[inArc].tailId;
	// 			sum += y[i][q][j];
	// 		}
	// 		obj += u_qj * sum * phi[q][j];
	// 	}
	// }
	// model.setObjective(obj, GRB_MINIMIZE);
	// model.update();
// }

Cut GuroSolver::solveSubProblem(const vector<vector<vector<shi>>> &y_bar) {
	// TODO: run all scenarios and return feasiblity/optimality cut.
	CutCoefficients coefficients;
	double rhs = 0;
	bool first = true;

	for (int i = 0; i < networkPtr->nScenarios; i++) {
		const auto cut = solveSubProblemInstance(y_bar);
		if (cut.cutType == FEASIBILITY) {
			// break and return feasibility cut.
			// cout << "Returning feasibility cut" << endl;
			return cut;
		}
		else {
			// aggregate this cut.
			if (first) {
				for (const auto& [k,v] : cut.cutCoeff) coefficients[k] = 0;
				first = false;
			}
			for (const auto& [k,v] : cut.cutCoeff) coefficients[k] += v;
			rhs += cut.RHS;
		}
		// cout << "Solved " <<i<< " th subproblem" << endl;
	}
	// cout << "Returning optimality cut" << endl;
	rhs /= networkPtr->nScenarios;
	for (auto& [_,v] : coefficients) v /= networkPtr->nScenarios;
	return {OPTIMALITY, rhs, coefficients};
}

double GuroSolver::solvepartialProblem(vector<int> partialSol, Network network) {
	auto n = network.n;
    double scNum = network.nScenarios * 1.0;

	cout << "network.nScenarios: " << network.nScenarios << endl;

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_IntParam_InfUnbdInfo, 1);

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
		   y_temp[j] = model.addVars(n, GRB_CONTINUOUS);
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

	for (int i = 0; i < n; i++) {
		for (int q = 0; q < n; q++) {
			for (int j = 0; j < n; j++) {
				model.addConstr(y[i][q][j] <= 1, "2cp");
			}
		}
	}
	model.update();

	for (size_t i = 0; i < partialSol.size(); i++) {
		auto decision = partialSol[i];
		if (decision != -1) {
			auto netArcId = network.processingOrder[i].second;
			uint iNetId = network.networkArcs[netArcId].tailId;
			uint qNetId = network.networkArcs[netArcId].headId;
			uint jNetId = network.networkArcs[decision].headId;
			model.addConstr(y[iNetId][qNetId][jNetId] == 1, "2dp");
		}else {
			auto netArcId = network.processingOrder[i].second;
			uint qNetId = network.networkArcs[netArcId].headId;
			uint iNetId = network.networkArcs[netArcId].tailId;
			for (auto outArc : network.networkNodes[qNetId].outgoingArcs) {
				uint jNetId = network.networkArcs[outArc].headId;
				model.addConstr(y[iNetId][qNetId][jNetId] == 0, "2dp");
			}
		}
	}

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

	cout << "zart o zurt" << endl;
	if(model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		cout << model.get(GRB_DoubleAttr_ObjVal) << endl;
		return model.get(GRB_DoubleAttr_ObjVal);
	}
	else {
		cout << "infeasable solution" << endl;
		return std::numeric_limits<int64_t>::min();
	}

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




}

Cut GuroSolver::solveSubProblemInstance(const vector<vector<vector<shi>>> &y_bar) {
	GRBModel model = GRBModel(environment);
	model.set(GRB_IntParam_Threads,1);
	model.set(GRB_IntParam_InfUnbdInfo, 1);

	const auto& nodes = networkPtr->networkNodes;
	const auto& arcs = networkPtr->networkArcs;
	double scNum = networkPtr->nScenarios;

	// var definition
	GRBVar* alpha = model.addVars(n, GRB_CONTINUOUS);
	GRBVar** beta = new GRBVar*[n];
	GRBVar** gamma = new GRBVar*[n];
	GRBVar** sigma = new GRBVar*[n];
	GRBVar** phi = new GRBVar*[n];
	GRBVar*** lambda = new GRBVar**[n];
	GRBVar*** mu = new GRBVar**[n];
	/// Initialilze variables ///
	for (int i = 0; i < n; i++)
	{
		beta[i] = model.addVars(n, GRB_CONTINUOUS);
		gamma[i] = model.addVars(n, GRB_CONTINUOUS);
		sigma[i] = model.addVars(n, GRB_CONTINUOUS);
		phi[i] = model.addVars(n, GRB_CONTINUOUS);
	}
	for (int i = 0; i < n; i++)
	{
		alpha[i].set(GRB_DoubleAttr_LB, -GRB_INFINITY);
		for (int q = 0; q < n; q++)
		{
			beta[i][q].set(GRB_DoubleAttr_LB, 0);
			gamma[i][q].set(GRB_DoubleAttr_LB, 0);
			sigma[i][q].set(GRB_DoubleAttr_LB, 0);
			phi[i][q].set(GRB_DoubleAttr_LB, 0);
		}
	}
	for (int i = 0; i < n; i++)
	{
		auto** lambda_temp = new GRBVar * [n];
		auto** mu_temp = new GRBVar * [n];
		for (int j = 0; j < n; j++)
		{
			lambda_temp[j] = model.addVars(n, GRB_CONTINUOUS);
			mu_temp[j] = model.addVars(n, GRB_CONTINUOUS);
		}
		lambda[i] = lambda_temp;
		mu[i] = mu_temp;
	}
	for (int i = 0; i < n; i++)
	{
		for (int q = 0; q < n; q++)
		{
			for (int j = 0; j < n; j++)
			{
				lambda[i][q][j].set(GRB_DoubleAttr_LB, 0);
				mu[i][q][j].set(GRB_DoubleAttr_LB, 0);
			}
		}
	}
	model.update();
	// model.optimize();



	/// Add constraints ///
	/// constraint number 1 & 2 ///
	for (uint arcID : networkPtr->A1) {
		GRBLinExpr LHS = 0;
		uint q = arcs[arcID].headId;
		uint i = arcs[arcID].tailId;
		int r_iq = arcs[arcID].rewards[0];
		LHS += alpha[q] - beta[i][q] + gamma[i][q];
		if (nodes[q].isVbar) {
			for (uint outArc : nodes[q].outgoingArcs) {
				uint j = arcs[outArc].headId;
				LHS += lambda[i][q][j] - mu[i][q][j];
			}
			LHS += sigma[i][q];
			model.addConstr(LHS >= r_iq, "7b");
		}else {
			model.addConstr(LHS >= r_iq, "7c");
		}
	}
	/// constraint number 3 & 4 ///
	for (uint arcID : networkPtr->A2) {
		GRBLinExpr LHS = 0;
		uint j = arcs[arcID].headId;
		uint q = arcs[arcID].tailId;
		int r_qj = arcs[arcID].rewards[0];
		LHS += -alpha[q] - beta[q][j] + gamma[q][j];
		if (nodes[q].isVbar) {
			for (uint outArc : nodes[q].incomingArcs) {
				uint i = arcs[outArc].tailId;
				LHS += -lambda[i][q][j] + mu[i][q][j];
			}
			LHS += phi[q][j];
			model.addConstr(LHS >= r_qj, "7b");
		}else {
			model.addConstr(LHS >= r_qj, "7e");
		}
	}
	//// constraint number 5 & 6 & 7 & 8 ////
	for (uint arcID : networkPtr->A3) {
		GRBLinExpr LHS = 0;
		uint j = arcs[arcID].headId;
		uint q = arcs[arcID].tailId;
		int r_qj = arcs[arcID].rewards[0];
		LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
		if (nodes[q].isVbar) {
			if (nodes[j].isVbar) {
				for (uint inArc : nodes[q].incomingArcs) {
					uint i = arcs[inArc].tailId;
					LHS += mu[i][q][j] - lambda[i][q][j];
				}
				for (uint outArc : nodes[j].outgoingArcs) {
					uint i = arcs[outArc].headId;
					LHS += lambda[q][j][i] - mu[q][j][i];
				}
				LHS += phi[q][j]+sigma[q][j];
				model.addConstr(LHS >= r_qj, "7f");
			}else {
				for (uint inArc : nodes[q].incomingArcs) {
					uint i = arcs[inArc].tailId;
					LHS += mu[i][q][j] - lambda[i][q][j];
				}
				LHS +=phi[q][j];
				model.addConstr(LHS >= r_qj, "7g");
			}
		}else {
			if (nodes[j].isVbar) {
				for (uint outArc : nodes[j].outgoingArcs) {
					uint i = arcs[outArc].headId;
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
	// for (auto arcID1 : networkPtr->A4) {
	// 	GRBLinExpr LHS = 0;
	// 	uint q = arcs[arcID1].headId;
	// 	uint i = arcs[arcID1].tailId;
	// 	uint r_iq = arcs[arcID1].rewards[scenario];
	// 	LHS += -beta[i][q] + gamma[i][q];
	// 	model.addConstr(LHS >= r_iq, "7j");
	// }
	model.addConstr(alpha[0] == 0);
	model.addConstr(alpha[n-1] == 0);


	model.update();




	const auto& processingOrder = networkPtr->processingOrder;
	CutCoefficients Y_bar_coef;
	double rhs = 0.0;
	CutType type = OPTIMALITY;
	for (uint i0 = 0; i0 < processingOrder.size(); i0++){
		uint arcID = processingOrder[i0].second;
		uint i = arcs[arcID].tailId;
		uint q = arcs[arcID].headId;
		for (uint outArc : nodes[q].outgoingArcs){
			uint j = arcs[outArc].headId;
			Y_bar_coef[make_tuple(i, q, j)] = 0.0;
		}
	}

	for (auto scenario = 0 ; scenario < networkPtr->nScenarios ; scenario++) {
		// cout << scenario << " - ";
		/// Set Objective function ///
		GRBLinExpr obj = 0;
		// first term
		for (int q = 0; q < n; q++){
			for (uint outArc : nodes[q].outgoingArcs){
				auto j = arcs[outArc].headId;
				int u_qj = arcs[outArc].upperCapacities[scenario];
				int l_qj = arcs[outArc].lowerCapacities[scenario];
				obj += u_qj*gamma[q][j] - l_qj*beta[q][j];
			}
		}
		// second term
		for (uint q : networkPtr->Vbar){
			for (uint inArc : nodes[q].incomingArcs){
				uint i = arcs[inArc].tailId;
				int u_iq = arcs[inArc].upperCapacities[scenario];
				for (uint outArc : nodes[q].outgoingArcs) {
					//uint i = arcs[inArc].tailId;
					uint j = arcs[outArc].headId;
					int u_qj = arcs[outArc].upperCapacities[scenario];
					obj += u_iq * (1-y_bar[i][q][j]) * lambda[i][q][j] + u_qj * (1-y_bar[i][q][j]) * mu[i][q][j]; // second term
					// sum += y_bar[i][q][j];
				}
				// third term here.
				// obj += u_iq * sum * sigma[i][q];
			}
		}
		// third term
		for (uint q : networkPtr->Vbar) {
			for (uint inArc : nodes[q].incomingArcs) {
				uint i = arcs[inArc].tailId;
				int u_iq = arcs[inArc].upperCapacities[scenario];
				int sum = 0;
				for (uint outArc : nodes[q].outgoingArcs) {
					uint j = arcs[outArc].headId;
					sum += y_bar[i][q][j];
				}
				obj += u_iq * sum * sigma[i][q];
			}
		}

		// fourth term
		for (uint q : networkPtr->Vbar){
			for (uint outArc : nodes[q].outgoingArcs){
				auto j = arcs[outArc].headId;
				int u_qj = arcs[outArc].upperCapacities[scenario];
				int sum = 0;
				for (auto inArc : nodes[q].incomingArcs){
					uint i = arcs[inArc].tailId;
					sum += y_bar[i][q][j];
				}
				obj += u_qj * sum * phi[q][j];
			}
		}
		model.setObjective(obj, GRB_MINIMIZE);
		model.update();
		model.optimize();

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL) {
			// first term
			for (uint q = 0; q < n ; q++){
				for (uint arcID : nodes[q].outgoingArcs){
					uint j = arcs[arcID].headId;
					rhs += (arcs[arcID].upperCapacities[scenario]/ scNum) * gamma[q][j].get(GRB_DoubleAttr_X)  ;
					rhs -= (arcs[arcID].lowerCapacities[scenario]/ scNum) * beta[q][j].get(GRB_DoubleAttr_X)  ;
				}
			}
			/// second term ///
			for (uint q : networkPtr->Vbar){
				for (uint inArcID : nodes[q].incomingArcs){
					uint i = arcs[inArcID].tailId;
					for (uint outArcID : nodes[q].outgoingArcs){
						uint j = arcs[outArcID].headId;
						rhs += (arcs[inArcID].upperCapacities[scenario] / scNum) * lambda[i][q][j].get(GRB_DoubleAttr_X) ;
						rhs += (arcs[outArcID].upperCapacities[scenario]/ scNum) * mu[i][q][j].get(GRB_DoubleAttr_X) ;
						Y_bar_coef[make_tuple(i, q, j)] -= (arcs[inArcID].upperCapacities[scenario] / scNum) * lambda[i][q][j].get(GRB_DoubleAttr_X) ;
						Y_bar_coef[make_tuple(i, q, j)] -= (arcs[outArcID].upperCapacities[scenario]/ scNum) * mu[i][q][j].get(GRB_DoubleAttr_X) ;
					}
				}
			}
			/// third term ///
			for (uint q : networkPtr->Vbar){
				for (uint inArcID : nodes[q].incomingArcs){
					uint i = arcs[inArcID].tailId;
					uint u_iq = arcs[inArcID].upperCapacities[scenario];
					for (uint outArcID : nodes[q].outgoingArcs){
						uint j = arcs[outArcID].headId;
						Y_bar_coef[make_tuple(i, q, j)] +=  (u_iq / scNum) * sigma[i][q].get(GRB_DoubleAttr_X) ;
					}
				}
			}
			/// fourth term ///
			for (uint q : networkPtr->Vbar){
				for (uint outArcID : nodes[q].outgoingArcs){
					uint j = arcs[outArcID].headId;
					uint u_qj = arcs[outArcID].upperCapacities[scenario];
					for (uint inArcID : nodes[q].incomingArcs){
						uint i = arcs[inArcID].tailId;
						Y_bar_coef[make_tuple(i, q, j)] += (u_qj  / scNum) * phi[q][j].get(GRB_DoubleAttr_X) ;
					}
				}
			}
			model.reset();
		}else {
			//////////////////////////////
			/////   feasibility cut   /////
			//////////////////////////////
			// CutCoefficients Y_bar_coef;
			double rhs = 0;
			for (uint i0 = 0; i0 < processingOrder.size(); i0++){
				uint arcID = processingOrder[i0].second;
				uint i = arcs[arcID].tailId;
				uint q = arcs[arcID].headId;
				for (uint outArc : nodes[q].outgoingArcs){
					uint j = arcs[outArc].headId;
					Y_bar_coef[make_tuple(i, q, j)] = 0;
				}
			}

			/// first term ///
			for (uint q = 0; q < n; q++){
				for (uint arcID : nodes[q].outgoingArcs){
					uint j = arcs[arcID].headId;
					rhs += arcs[arcID].upperCapacities[scenario] * gamma[q][j].get(GRB_DoubleAttr_UnbdRay);
					rhs -= arcs[arcID].lowerCapacities[scenario] * beta[q][j].get(GRB_DoubleAttr_UnbdRay);
				}
			}

			/// second term ///
			for (uint q : networkPtr->Vbar){
				for (uint inArcID : nodes[q].incomingArcs){
					uint i = arcs[inArcID].tailId;
					for (uint outArcID : nodes[q].outgoingArcs){
						uint j = arcs[outArcID].headId;
						rhs += arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
						rhs += arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
						Y_bar_coef[make_tuple(i, q, j)] -= arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
						Y_bar_coef[make_tuple(i, q, j)] -= arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					}
				}
			}

			/// third term ///
			for (uint q : networkPtr->Vbar){
				for (uint inArcID : nodes[q].incomingArcs){
					uint i = arcs[inArcID].tailId;
					uint u_iq = arcs[inArcID].upperCapacities[scenario];
					for (uint outArcID : nodes[q].outgoingArcs){
						uint j = arcs[outArcID].headId;
						Y_bar_coef[make_tuple(i, q, j)] += u_iq * sigma[i][q].get(GRB_DoubleAttr_UnbdRay);
					}
				}
			}

			/// fourth term ///
			for (uint q : networkPtr->Vbar){
				for (uint outArcID : nodes[q].outgoingArcs){
					uint j = arcs[outArcID].headId;
					uint u_qj = arcs[outArcID].upperCapacities[scenario];
					for (uint inArcID : nodes[q].incomingArcs){
						uint i = arcs[inArcID].tailId;
						Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_UnbdRay);
					}
				}
			}
			cout << "-f-" ;
			type = FEASIBILITY;
			return Cut{type, rhs, Y_bar_coef};
		}
	}
	// std::setprecision(15);
	cout <<  std::setprecision(15) << "";
	cout << "-o-" ;
	return Cut{type, rhs, Y_bar_coef};
}

















	// if (status == GRB_OPTIMAL){
		//////////////////////////
		/// optimality cut ///////
		/////////////////////////
		//CutCoefficients Y_bar_coef;
		// for (uint i0 = 0; i0 < processingOrder.size(); i0++){
		// 	uint arcID = processingOrder[i0].second;
		// 	uint i = arcs[arcID].tailId;
		// 	uint q = arcs[arcID].headId;
		// 	for (uint outArc : nodes[q].outgoingArcs){
		// 		uint j = arcs[outArc].headId;
		// 		Y_bar_coef[make_tuple(i, q, j)] = 0;
		// 	}
		// }
		// // double rhs = 0;

		/// first term ///
		// for (uint q = 0; q < n; q++){
		// 	for (uint arcID : nodes[q].outgoingArcs){
		// 		uint j = arcs[arcID].headId;
		// 		rhs += arcs[arcID].upperCapacities[scenario] * gamma[q][j].get(GRB_DoubleAttr_X);
		// 		rhs -= arcs[arcID].lowerCapacities[scenario] * beta[q][j].get(GRB_DoubleAttr_X);
		// 	}
		// }

		// /// second term ///
		// for (uint q : networkPtr->Vbar){
		// 	for (uint inArcID : nodes[q].incomingArcs){
		// 		uint i = arcs[inArcID].tailId;
		// 		for (uint outArcID : nodes[q].outgoingArcs){
		// 			uint j = arcs[outArcID].headId;
		// 			rhs += arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_X);
		// 			rhs += arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_X);
		// 			Y_bar_coef[make_tuple(i, q, j)] -= arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_X);
		// 			Y_bar_coef[make_tuple(i, q, j)] -= arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_X);
		// 		}
		// 	}
		// }

		// /// third term ///
		// for (uint q : networkPtr->Vbar){
		// 	for (uint inArcID : nodes[q].incomingArcs){
		// 		uint i = arcs[inArcID].tailId;
		// 		uint u_iq = arcs[inArcID].upperCapacities[scenario];
		// 		for (uint outArcID : nodes[q].outgoingArcs){
		// 			uint j = arcs[outArcID].headId;
		// 			Y_bar_coef[make_tuple(i, q, j)] += u_iq * sigma[i][q].get(GRB_DoubleAttr_X);
		// 		}
		// 	}
		// }

		// /// fourth term ///
		// for (uint q : networkPtr->Vbar){
		// 	for (uint outArcID : nodes[q].outgoingArcs){
		// 		uint j = arcs[outArcID].headId;
		// 		uint u_qj = arcs[outArcID].upperCapacities[scenario];
		// 		for (uint inArcID : nodes[q].incomingArcs){
		// 			uint i = arcs[inArcID].tailId;
		// 			Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_X);
		// 		}
		// 	}
		// }
		// model.reset(1);
		// return Cut{OPTIMALITY, rhs, Y_bar_coef};
		// type = OPTIMALITY;
// 	}
// 	else {
// 		//////////////////////////////
// 		/////   feasibility cut   /////
// 		//////////////////////////////
//
// 		// CutCoefficients Y_bar_coef;
// 		for (uint i0 = 0; i0 < processingOrder.size(); i0++){
// 			uint arcID = processingOrder[i0].second;
// 			uint i = arcs[arcID].tailId;
// 			uint q = arcs[arcID].headId;
// 			for (uint outArc : nodes[q].outgoingArcs){
// 				uint j = arcs[outArc].headId;
// 				// cout << "(" << i << "," << q << "," << j << ")" << endl;
// 				Y_bar_coef[make_tuple(i, q, j)] = 0;
// 			}
// 		}
// 		// double rhs = 0;
//
// 		/// first term ///
// 		for (uint q = 0; q < n; q++){
// 			for (uint arcID : nodes[q].outgoingArcs){
// 				uint j = arcs[arcID].headId;
// 				rhs += arcs[arcID].upperCapacities[scenario] * gamma[q][j].get(GRB_DoubleAttr_UnbdRay);
// 				rhs -= arcs[arcID].lowerCapacities[scenario] * beta[q][j].get(GRB_DoubleAttr_UnbdRay);
// 			}
// 		}
//
// 		/// second term ///
// 		for (uint q : networkPtr->Vbar){
// 			for (uint inArcID : nodes[q].incomingArcs){
// 				uint i = arcs[inArcID].tailId;
// 				for (uint outArcID : nodes[q].outgoingArcs){
// 					uint j = arcs[outArcID].headId;
// 					rhs += arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
// 					rhs += arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
// 					Y_bar_coef[make_tuple(i, q, j)] -= arcs[inArcID].upperCapacities[scenario] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
// 					Y_bar_coef[make_tuple(i, q, j)] -= arcs[outArcID].upperCapacities[scenario] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
// 				}
// 			}
// 		}
//
// 		/// third term ///
// 		for (uint q : networkPtr->Vbar){
// 			for (uint inArcID : nodes[q].incomingArcs){
// 				uint i = arcs[inArcID].tailId;
// 				uint u_iq = arcs[inArcID].upperCapacities[scenario];
// 				for (uint outArcID : nodes[q].outgoingArcs){
// 					uint j = arcs[outArcID].headId;
// 					Y_bar_coef[make_tuple(i, q, j)] += u_iq * sigma[i][q].get(GRB_DoubleAttr_UnbdRay);
// 				}
// 			}
// 		}
//
// 		/// fourth term ///
// 		for (uint q : networkPtr->Vbar){
// 			for (uint outArcID : nodes[q].outgoingArcs){
// 				uint j = arcs[outArcID].headId;
// 				uint u_qj = arcs[outArcID].upperCapacities[scenario];
// 				for (uint inArcID : nodes[q].incomingArcs){
// 					uint i = arcs[inArcID].tailId;
// 					Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_UnbdRay);
// 				}
// 			}
// 		}
// 		// model.reset(1);
// 		// return Cut{FEASIBILITY, rhs, Y_bar_coef};
// 		type = FEASIBILITY;
// 	}
// 	// reset the model.
// 	//model.reset(1); // clean up gurobi variables.
// 	for (int i = 0; i < n; i++){
// 		delete[](beta[i]);
// 		delete[](gamma[i]);
// 		delete[](sigma[i]);
// 		delete[](phi[i]);
//
// 		for (int j = 0; j < n; j++){
// 			delete[](lambda[i][j]);
// 			delete[](mu[i][j]);
// 		}
// 		delete[](lambda[i]);
// 		delete[](mu[i]);
// 	}
//
// 	delete[](alpha);
// 	delete[](beta);
// 	delete[](gamma);
// 	delete[](sigma);
// 	delete[](phi);
// 	delete[](lambda);
// 	delete[](mu);
// 	// return cut,
// 	return {type, rhs, Y_bar_coef};
// }