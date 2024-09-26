//
// Created by erfank on 9/24/2024.
//

#include "CGO.h"
#include "C:\gurobi911\win64\include\gurobi_c++.h"
#include <iostream>
#include <ostream>
#include <map>
using namespace std;
void generateCut(vector<int> &W_solution, Network &network , Cut &newCut , int s) {
    cout << "Cut Generation Started" << endl;
    ////////////////////////////
    /// from w solution to y ///
    ////////////////////////////
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
	cout << "Solution transformation finished" << endl;
	for (int i0 = 0; i0 < network.n; i0++) {
		for (int i1 = 0; i1 < network.n; i1++) {
			for (int i2 = 0; i2 < network.n; i2++) {
				if (y_bar[i0][i1][i2] == 1) {
					cout << i0 << " - " << i1 << " - " << i2 << endl;
				}
			}
		}
	}
    ///////////////////////////////////////////
    /// Gurobi env creation & Var definition///
    ///////////////////////////////////////////
    GRBEnv env = GRBEnv();
    GRBModel Dual_subproblem = GRBModel(env);
    Dual_subproblem.set(GRB_IntParam_InfUnbdInfo, 1);

    //GRBEnv environment;
    // vector<GRBVar> alpha;
    // vector<vector<GRBVar>> beta;
    // vector<vector<GRBVar>> gamma;
    // vector<vector<vector<GRBVar>>> lambda;
    // vector<vector<vector<GRBVar>>> mu;
    // vector<vector<GRBVar>> sigma;
    // vector<vector<GRBVar>> phi;
	cout << "variable definition finished " << endl;
    uint n = network.n;

  //   for (int i = 0; i < n; i++) {
  //       alpha[i] = Dual_subproblem.addVar(-GRB_INFINITY, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //
  //       beta[i] = vector<GRBVar>(n);beta[i] = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
  //       gamma[i] = vector<GRBVar>(n);
  //       sigma[i] = vector<GRBVar>(n);
  //       phi[i] = vector<GRBVar>(n);
  //
  //       lambda[i] = vector<vector<GRBVar>>(n);
  //       mu[i] = vector<vector<GRBVar>>(n);
		// cout << "interation number: " << i << endl;
  //       for (int j = 0; j < n; j++) {
  //           beta[i][j] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //           gamma[i][j] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //           sigma[i][j] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //           phi[i][j] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //
  //           lambda[i][j] = vector<GRBVar>(n);
  //           mu[i][j] = vector<GRBVar>(n);
  //
  //           for (int k = 0; k < n; k++) {
  //               lambda[i][j][k] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //               mu[i][j][k] = Dual_subproblem.addVar(0, GRB_INFINITY, NULL, GRB_CONTINUOUS);
  //           }
  //       }
  //   }
	GRBVar** beta = new GRBVar * [n];
	GRBVar** gamma = new GRBVar * [n];
	GRBVar*** lambda = new GRBVar * *[n];
	GRBVar*** mu = new GRBVar * *[n];
	GRBVar** sigma = new GRBVar * [n];
	GRBVar** phi = new GRBVar * [n];
	GRBVar* alpha = Dual_subproblem.addVars(n, GRB_CONTINUOUS);
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
		GRBVar** lambda_temp = new GRBVar * [n];
		GRBVar** mu_temp = new GRBVar * [n];
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
    for (int q = 0; q < n; q++) {
        for (uint outArc : network.networkNodes[q].outgoingArcs)
        {
            uint j = network.networkArcs[outArc].headId;
            Dual_objective -= network.networkArcs[outArc].lowerCapacities[s] * beta[q][j];
            Dual_objective += network.networkArcs[outArc].upperCapacities[s] * gamma[q][j];
        }
    }
    /// second term ///
    for (uint q: network.Vbar)
    {
        for (uint inArc : network.networkNodes[q].incomingArcs)
        {
            for (uint outArc : network.networkNodes[q].outgoingArcs)
            {
                uint i = network.networkArcs[inArc].tailId;
                uint j = network.networkArcs[outArc].headId;
                Dual_objective += network.networkArcs[inArc].upperCapacities[s] * (1 - y_bar[i][q][j]) * lambda[i][q][j];
                Dual_objective += network.networkArcs[outArc].upperCapacities[s] * (1 - y_bar[i][q][j]) * mu[i][q][j];

            }
        }
    }
    /// third term ///
    for (uint q: network.Vbar)
    {
        for (uint inArc : network.networkNodes[q].incomingArcs)
        {
            uint i = network.networkArcs[inArc].tailId;
            uint u_iq = network.networkArcs[inArc].upperCapacities[s];
            for (uint outArc : network.networkNodes[q].outgoingArcs)
            {
                uint j = network.networkArcs[outArc].headId;
                Dual_objective += u_iq * y_bar[i][q][j] * sigma[i][q];
            }
        }
    }
    /// fourth term ///
    for (uint q : network.Vbar)
    {
        for (uint outArc : network.networkNodes[q].outgoingArcs)
        {
            uint j = network.networkArcs[outArc].headId;
            uint u_qj = network.networkArcs[outArc].upperCapacities[s];
            for (uint inArc : network.networkNodes[q].incomingArcs)
            {
                uint i = network.networkArcs[inArc].tailId;
                Dual_objective += u_qj * y_bar[i][q][j] * phi[q][j];
            }
        }
    }
    Dual_subproblem.setObjective(Dual_objective, GRB_MINIMIZE);

    ///////////////////////////////
    ///       constraints       ///
    ///////////////////////////////
    ///
    /// constraint number 1 & 2 ///
    for (uint arcID1 : network.A1)
    {
        GRBLinExpr LHS = 0;
        uint i = network.networkArcs[arcID1].tailId;
        uint q = network.networkArcs[arcID1].headId;
        uint r_iq = network.networkArcs[arcID1].rewards[s];
        if (network.isNodeInVbar[q])
        {
            LHS += alpha[q] - beta[i][q] + gamma[i][q];
            for (uint arcID2 : network.networkNodes[q].outgoingArcs)
            {
                uint j = network.networkArcs[arcID2].headId;
                LHS += lambda[i][q][j] - mu[i][q][j];
            }
            LHS += sigma[i][q];
            Dual_subproblem.addConstr(LHS >= r_iq, "7b");
        }
        else
        {
            LHS += alpha[q] - beta[i][q] + gamma[i][q];
            Dual_subproblem.addConstr(LHS >= r_iq, "7c");
        }
    }

    /// constraint number 3 & 4 ///

    for (uint arcID1 : network.A2)
    {
        GRBLinExpr LHS = 0;
        uint j = network.networkArcs[arcID1].headId;
        uint q = network.networkArcs[arcID1].tailId;
        uint r_qj = network.networkArcs[arcID1].rewards[s];
        if (network.isNodeInVbar[q])
        {
            LHS += -alpha[q] - beta[q][j] + gamma[q][j];
            for (uint arcID2 : network.networkNodes[q].incomingArcs)
            {
                uint i = network.networkArcs[arcID2].tailId;
                LHS += -lambda[i][q][j] + mu[i][q][j];
            }
            LHS += phi[q][j];
            Dual_subproblem.addConstr(LHS >= r_qj, "7d");
        }
        else
        {
            LHS += -alpha[q] - beta[q][j] + gamma[q][j];
            Dual_subproblem.addConstr(LHS >= r_qj, "7e");
        }
    }

    //// constraint number 5 & 6 & 7 & 8 ////

    for (uint arcID1 : network.A3)
    {
        GRBLinExpr LHS = 0;
        uint j = network.networkArcs[arcID1].headId;
        uint q = network.networkArcs[arcID1].tailId;
        uint r_qj = network.networkArcs[arcID1].rewards[s];
        if (network.isNodeInVbar[q]) {
            if (network.isNodeInVbar[j]) {
                LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
                for (uint arcID2 : network.networkNodes[q].incomingArcs)
                {
                    uint i = network.networkArcs[arcID2].tailId;
                    LHS += -lambda[i][q][j] + mu[i][q][j];
                }

                for (uint arcID2 : network.networkNodes[j].outgoingArcs)
                {
                    uint i = network.networkArcs[arcID2].headId;
                    LHS += lambda[q][j][i] - mu[q][j][i];
                }
                LHS += sigma[q][j] + phi[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7f");
            }
            else {
                LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
                for (uint arcID2 : network.networkNodes[q].incomingArcs)
                {
                    uint i = network.networkArcs[arcID2].tailId;
                    LHS += -lambda[i][q][j] + mu[i][q][j];
                }
                LHS += phi[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7g");
            }
        }
        else
        {
            if (network.isNodeInVbar[j]) {
                LHS += -alpha[q] + alpha[j] + beta[q][j] + gamma[q][j];
                for (uint arcID2 : network.networkNodes[j].outgoingArcs)
                {
                    uint i = network.networkArcs[arcID2].headId;
                    LHS += lambda[q][j][i] - mu[q][j][i];
                }
                LHS += sigma[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7h");
            }
            else {
                LHS += -alpha[q] + alpha[j] - beta[q][j] + gamma[q][j];
                Dual_subproblem.addConstr(LHS >= r_qj, "7i");
            }
        }
    }

    //// constraint number 9 ////

    for (int arcID1 : network.A4)
    {
        GRBLinExpr LHS = 0;
        uint q = network.networkArcs[arcID1].headId;
        uint i = network.networkArcs[arcID1].tailId;
        uint r_iq = network.networkArcs[arcID1].rewards[s];
        if (network.isNodeInVbar[q] == 0 ) {
            if (network.isNodeInVbar[i] == 0 ) {
                LHS += -beta[i][q] + gamma[i][q];
                Dual_subproblem.addConstr(LHS >= r_iq, "7i");
            }
        }
    }

    Dual_subproblem.update();
    Dual_subproblem.optimize();
    int dual_status = Dual_subproblem.get(GRB_IntAttr_Status);

    /////////////////////////////
    /////  CUT COMPUTATION  /////
    /////////////////////////////
    if(dual_status == GRB_OPTIMAL) {
        //////////////////////////////
        /////   optimality cut   /////
        //////////////////////////////

        map<tuple<int, int, int>, double> Y_bar_coef;
	    for (uint i0 = 0; i0 < network.processingOrder.size(); i0++)
	    {
		    uint arcID = network.processingOrder[i0].second;
		    uint i = network.networkArcs[arcID].tailId;
		    uint q = network.networkArcs[arcID].headId;
		    for (uint outArc : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[outArc].headId;
			    Y_bar_coef[make_tuple(i, q, j)] = 0;
		    }
	    }
	    double rhs = 0;

	    /// first term ///

	    for (uint q = 0; q < network.n; q++)
	    {
		    for (uint arcID : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[arcID].headId;
			    rhs += static_cast<double>(network.networkArcs[arcID].upperCapacities[s]) * gamma[q][j].get(GRB_DoubleAttr_X) - static_cast<double>(network.networkArcs[arcID].lowerCapacities[s]) * beta[q][j].get(GRB_DoubleAttr_X);
		    }
	    }

	    /// second term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint inArcID : network.networkNodes[q].incomingArcs)
		    {
			    uint i = network.networkArcs[inArcID].tailId;
			    for (uint outArcID : network.networkNodes[q].outgoingArcs)
			    {
				    uint j = network.networkArcs[outArcID].headId;
				    rhs += network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_X);
				    rhs += network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_X);
				    Y_bar_coef[make_tuple(i, q, j)] -= static_cast<double>(network.networkArcs[inArcID].upperCapacities[s]) * lambda[i][q][j].get(GRB_DoubleAttr_X);
				    Y_bar_coef[make_tuple(i, q, j)] -= static_cast<double>(network.networkArcs[outArcID].upperCapacities[s]) * mu[i][q][j].get(GRB_DoubleAttr_X);
			    }
		    }
	    }

	    /// third term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint inArcID : network.networkNodes[q].incomingArcs)
		    {
			    uint i = network.networkArcs[inArcID].tailId;
			    uint u_i_q = network.networkArcs[inArcID].upperCapacities[s];
			    for (uint outArcID : network.networkNodes[q].outgoingArcs)
			    {
				    uint j = network.networkArcs[outArcID].headId;
				    Y_bar_coef[make_tuple(i, q, j)] += static_cast<double>(u_i_q) * sigma[i][q].get(GRB_DoubleAttr_X);
			    }
		    }
	    }

	    /// fourth term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint outArcID : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[outArcID].headId;
			    uint u_q_j = network.networkArcs[outArcID].upperCapacities[s];
			    for (uint inArcID : network.networkNodes[q].incomingArcs)
			    {
				    uint i = network.networkArcs[inArcID].tailId;
				    Y_bar_coef[make_tuple(i, q, j)] += static_cast<double>(u_q_j) * phi[i][q].get(GRB_DoubleAttr_X);
			    }
		    }
	    }
    	newCut.type = 'o';
    	newCut.cutCoef = Y_bar_coef;
    	newCut.RHS = rhs;

    }else {

        //////////////////////////////
        /////   feasibilty cut   /////
        //////////////////////////////

        map<tuple<int, int, int>, double> Y_bar_coef;
	    for (uint i0 = 0; i0 < network.processingOrder.size(); i0++)
	    {
		    uint arcID = network.processingOrder[i0].second;
		    uint i = network.networkArcs[arcID].tailId;
		    uint q = network.networkArcs[arcID].headId;
		    for (uint outArc : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[outArc].headId;
			    Y_bar_coef[make_tuple(i, q, j)] = 0;
		    }
	    }
	    double rhs = 0;

	    /// first term ///

	    for (uint q = 0; q < network.n; q++)
	    {
		    for (uint arcID : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[arcID].headId;
			    rhs += static_cast<double>(network.networkArcs[arcID].upperCapacities[s]) * gamma[q][j].get(GRB_DoubleAttr_UnbdRay) - static_cast<double>(network.networkArcs[arcID].lowerCapacities[s]) * beta[q][j].get(GRB_DoubleAttr_UnbdRay);
		    }
	    }

	    /// second term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint inArcID : network.networkNodes[q].incomingArcs)
		    {
			    uint i = network.networkArcs[inArcID].tailId;
			    for (uint outArcID : network.networkNodes[q].outgoingArcs)
			    {
				    uint j = network.networkArcs[outArcID].headId;
				    rhs += network.networkArcs[inArcID].upperCapacities[s] * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
				    rhs += network.networkArcs[outArcID].upperCapacities[s] * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
				    Y_bar_coef[make_tuple(i, q, j)] -= static_cast<double>(network.networkArcs[inArcID].upperCapacities[s]) * lambda[i][q][j].get(GRB_DoubleAttr_UnbdRay);
				    Y_bar_coef[make_tuple(i, q, j)] -= static_cast<double>(network.networkArcs[outArcID].upperCapacities[s]) * mu[i][q][j].get(GRB_DoubleAttr_UnbdRay);
			    }
		    }
	    }

	    /// third term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint inArcID : network.networkNodes[q].incomingArcs)
		    {
			    uint i = network.networkArcs[inArcID].tailId;
			    uint u_i_q = network.networkArcs[inArcID].upperCapacities[s];
			    for (uint outArcID : network.networkNodes[q].outgoingArcs)
			    {
				    uint j = network.networkArcs[outArcID].headId;
				    Y_bar_coef[make_tuple(i, q, j)] += static_cast<double>(u_i_q) * sigma[i][q].get(GRB_DoubleAttr_UnbdRay);
			    }
		    }
	    }

	    /// fourth term ///

	    for (uint q : network.Vbar)
	    {
		    for (uint outArcID : network.networkNodes[q].outgoingArcs)
		    {
			    uint j = network.networkArcs[outArcID].headId;
			    uint u_q_j = network.networkArcs[outArcID].upperCapacities[s];
			    for (uint inArcID : network.networkNodes[q].incomingArcs)
			    {
				    uint i = network.networkArcs[inArcID].tailId;
				    Y_bar_coef[make_tuple(i, q, j)] += static_cast<double>(u_q_j) * phi[i][q].get(GRB_DoubleAttr_UnbdRay);
			    }
		    }
	    }

    	newCut.type = 'f';
    	newCut.cutCoef = Y_bar_coef;
    	newCut.RHS = rhs;
    }






















}