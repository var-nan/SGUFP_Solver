//
// Created by nandgate on 10/25/2024.
//

#ifndef NEWCGO_H
#define NEWCGO_H
#include "gurobi_c++.h"
#include <iostream>
#include "Cut.h"
#include "Network.h"

inline Cut generateNewCut(const Network& network, vector<int>& W_solution) {
	vector<int> y_1(network.n, 0);
	vector<vector<int>> y_2(network.n, y_1);
	vector < vector<vector<int>>> y_bar(network.n, y_2);
	for (int i0 = 0; i0 < W_solution.size(); i0++)
	{
		if (W_solution[i0] != -1){
			uint arcID = network.processingOrder[i0].second;
			uint q = network.networkArcs[arcID].headId;
			uint i = network.networkArcs[arcID].tailId;
			uint j = network.networkArcs[W_solution[i0]].headId;
			y_bar[i][q][j] = 1;
		}
	}
	int s = 0;
    GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 0);
    GRBModel Dual_subproblem = GRBModel(env);
    Dual_subproblem.set(GRB_IntParam_InfUnbdInfo, 1);

    cout << "variable definition finished " << endl;
    int n = static_cast<int>(network.n);
    auto* alpha = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
    auto** beta = new GRBVar * [n];
    auto** gamma = new GRBVar * [n];
    auto** sigma = new GRBVar * [n];
    auto** phi = new GRBVar * [n];
    auto*** lambda = new GRBVar * *[n];
    auto*** mu = new GRBVar * *[n];
    for (int i = 0; i < n; i++)
    {
        beta[i] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
        gamma[i] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
        sigma[i] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
        phi[i] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
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
            lambda_temp[j] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
            mu_temp[j] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
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

    cout << "variables populated.. " << endl;
    Dual_subproblem.update();



    //////////////////////////
    /// objective function ///
    //////////////////////////

    GRBLinExpr Dual_objective = 0;
    /// first term ///
    for (uint q = 0 ; q < n ; q++){
      for (uint outArc : network.networkNodes[q].outgoingArcs) {
          uint j = network.networkArcs[outArc].headId;
          int u_qj = network.networkArcs[outArc].upperCapacities[s];
          int l_qj = network.networkArcs[outArc].lowerCapacities[s];
          Dual_objective += u_qj * gamma[q][j];
          Dual_objective -= l_qj * beta[q][j];
      }
    }
    /// second term ///
    for (uint q : network.Vbar) {
        for (uint inArc : network.networkNodes[q].incomingArcs) {
            for (uint outArc : network.networkNodes[q].outgoingArcs) {
                uint j = network.networkArcs[outArc].headId;
                uint i = network.networkArcs[inArc].tailId;
                int u_qj = network.networkArcs[outArc].upperCapacities[s];
                int u_iq = network.networkArcs[inArc].upperCapacities[s];
                Dual_objective += u_iq * (1-y_bar[i][q][j]) * lambda[i][q][j];
                // Dual_objective -= u_iq * y_bar[i][q][j] * lambda[i][q][j];
                Dual_objective += u_qj * (1-y_bar[i][q][j]) * mu[i][q][j];
                // Dual_objective -= u_qj * y_bar[i][q][j] * mu[i][q][j];
            }
        }
    }
    /// third term ///
    for (uint q : network.Vbar) {
        for (uint inArc : network.networkNodes[q].incomingArcs) {
            uint i = network.networkArcs[inArc].tailId;
            int u_iq= network.networkArcs[inArc].upperCapacities[s];
            int sum = 0;
            for (uint outArc : network.networkNodes[q].outgoingArcs) {
                uint j = network.networkArcs[outArc].headId;
                sum += y_bar[i][q][j];
            }
            Dual_objective += u_iq * sum * sigma[i][q];
        }
    }
    /// fourth term ///
    for (uint q : network.Vbar) {
        for (uint outArc : network.networkNodes[q].outgoingArcs) {
            uint j = network.networkArcs[outArc].headId;
            int u_qj= network.networkArcs[outArc].upperCapacities[s];
            int sum = 0;
            for (uint inArc : network.networkNodes[q].incomingArcs) {
                uint i = network.networkArcs[inArc].tailId;
                sum += y_bar[i][q][j];
            }
            Dual_objective += u_qj * sum * phi[q][j];
        }
    }

    Dual_subproblem.setObjective(Dual_objective, GRB_MINIMIZE);
    Dual_subproblem.update();

    ///////////////////////////////
    ///       constraints       ///
    ///////////////////////////////

    /// constraint number 1 & 2 ///
    for (uint arcID : network.A1) {
        GRBLinExpr LHS = 0;
        uint q = network.networkArcs[arcID].headId;
        uint i = network.networkArcs[arcID].tailId;
        int r_iq = network.networkArcs[arcID].rewards[s];
        LHS += alpha[q] - beta[i][q] + gamma[i][q];
        if (network.networkNodes[q].isVbar) {
            for (uint outArc : network.networkNodes[q].outgoingArcs) {
                uint j = network.networkArcs[outArc].headId;
                LHS += lambda[i][q][j] - mu[i][q][j];
            }
            LHS += sigma[i][q];
            Dual_subproblem.addConstr(LHS >= r_iq, "7b");
        }else {
            Dual_subproblem.addConstr(LHS >= r_iq, "7c");
        }
    }

    /// constraint number 3 & 4 ///
    for (uint arcID : network.A2) {
        GRBLinExpr LHS = 0;
        uint j = network.networkArcs[arcID].headId;
        uint q = network.networkArcs[arcID].tailId;
        int r_qj = network.networkArcs[arcID].rewards[s];
        LHS += -alpha[q] - beta[q][j] + gamma[q][j];
        if (network.networkNodes[q].isVbar) {
            for (uint incArc : network.networkNodes[q].incomingArcs) {
                uint i = network.networkArcs[incArc].tailId;
                LHS += -lambda[i][q][j] + mu[i][q][j];
            }
            LHS += phi[q][j];
            Dual_subproblem.addConstr(LHS >= r_qj, "7b");
        }else {
            Dual_subproblem.addConstr(LHS >= r_qj, "7e");
        }
    }

    //// constraint number 5 & 6 & 7 & 8 ////

    for (uint arcID : network.A3) {
        GRBLinExpr LHS = 0;
        uint j = network.networkArcs[arcID].headId;
        uint q = network.networkArcs[arcID].tailId;
        int r_qj = network.networkArcs[arcID].rewards[s];
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
                Dual_subproblem.addConstr(LHS >= r_qj, "7f");
            }else {
                for (uint inArc : network.networkNodes[q].incomingArcs) {
                    uint i = network.networkArcs[inArc].tailId;
                    LHS += mu[i][q][j] - lambda[i][q][j];
                }
                LHS +=phi[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7g");
            }
        }else {
            if (network.networkNodes[j].isVbar) {
                for (uint outArc : network.networkNodes[j].outgoingArcs) {
                    uint i = network.networkArcs[outArc].headId;
                    LHS += lambda[q][j][i] - mu[q][j][i];
                }
                LHS += sigma[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7h");
            }else {
                Dual_subproblem.addConstr(LHS >= r_qj, "7i");
            }
        }
    }

    //// constraint number 9 ////

    // for (int arcID1 : network.A4)
    // {
    //     GRBLinExpr LHS = 0;
    //     uint q = network.networkArcs[arcID1].headId;
    //     uint i = network.networkArcs[arcID1].tailId;
    //     uint r_iq = network.networkArcs[arcID1].rewards[s];
    //     LHS += -beta[i][q] + gamma[i][q];
    //     Dual_subproblem.addConstr(LHS >= r_iq, "7j");
    // }
    //
    // for (int q=0 ; q<n ; q++) {
    //     if (q== 0 || q==n-1) {
    //         Dual_subproblem.addConstr(alpha[q] == 0, "7k");
    //     }
    // }
    Dual_subproblem.update();
    Dual_subproblem.optimize();
    int dual_status = Dual_subproblem.get(GRB_IntAttr_Status);

	if(dual_status == GRB_OPTIMAL) {
		//////////////////////////////
		/////   optimality cut   /////
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
				rhs += network.networkArcs[arcID].upperCapacities[s] * gamma[q][j].get(GRB_DoubleAttr_X);
				rhs -= network.networkArcs[arcID].lowerCapacities[s] * beta[q][j].get(GRB_DoubleAttr_X);
			}
		}

		/// second term ///

		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					rhs += network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_X);
					rhs += network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_X);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_X);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_X);
				}
			}
		}

		/// third term ///

		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
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
				uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
				for (uint inArcID : network.networkNodes[q].incomingArcs){
					uint i = network.networkArcs[inArcID].tailId;
					Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_X);
				}
			}
		}

		Cut newCut{OPTIMALITY, rhs, Y_bar_coef};

		// newCut.type = 'o';
		// newCut.cutCoef = Y_bar_coef;
		// newCut.RHS = rhs;
		return newCut;
	}
	else {

		//////////////////////////////
		/////   feasibilty cut   /////
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
				rhs += network.networkArcs[arcID].upperCapacities[s] * gamma[q][j].get(GRB_DoubleAttr_UnbdRay);
				rhs -= network.networkArcs[arcID].lowerCapacities[s] * beta[q][j].get(GRB_DoubleAttr_UnbdRay);
			}
		}

		/// second term ///

		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				for (uint outArcID : network.networkNodes[q].outgoingArcs){
					uint j = network.networkArcs[outArcID].headId;
					rhs += network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					rhs += network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
					Y_bar_coef[make_tuple(i, q, j)] -= network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
				}
			}
		}

		/// third term ///

		for (uint q : network.Vbar){
			for (uint inArcID : network.networkNodes[q].incomingArcs){
				uint i = network.networkArcs[inArcID].tailId;
				uint u_iq = network.networkArcs[inArcID].upperCapacities[s];
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
				uint u_qj = network.networkArcs[outArcID].upperCapacities[s];
				for (uint inArcID : network.networkNodes[q].incomingArcs){
					uint i = network.networkArcs[inArcID].tailId;
					Y_bar_coef[make_tuple(i, q, j)] += u_qj * phi[q][j].get(GRB_DoubleAttr_UnbdRay);
				}
			}
		}

		// vector<double> lowerBounds(network.processingOrder.size());
		// const auto& coeff = Y_bar_coef;
		//
		// for (const auto[id_t, arcId]: network.processingOrder){
		// 	const auto& arc = network.networkArcs[arcId];
		// 	auto i = arc.tailId;
		// 	auto q = arc.headId;
		//
		// 	const auto& node = network.networkNodes[q];
		// 	int max = 0;
		// 	for (const auto j: node.outNodeIds){
		// 		auto c = coeff.at(make_tuple(static_cast<int>(i),static_cast<int>(q),static_cast<int>(j)));
		// 		if (c > max) max = c;
		// 	}
		// 	lowerBounds[id_t] = max;
		// }
		// int pref = 0;
		// for (auto start = lowerBounds.rbegin(); start != lowerBounds.rend(); start++){
		// 	*start = *start + pref;
		// 	pref = *start;
		// }
		// lowerBounds.push_back(0);

		Cut newCut{FEASIBILITY, rhs, Y_bar_coef};
		return newCut;

		// newCut.type = 'f';
		// newCut.cutCoef = Y_bar_coef;
		// newCut.RHS = rhs;
		// newCut.heuristic = lowerBounds;
	}
}
#endif //NEWCGO_H
