//
// Created by nandgate on 10/27/2024.
//

#ifndef NODEEXPLORER_H
#define NODEEXPLORER_H

#include "DD.h"
#include "Cut.h"
#include "gurobi_c++.h"
#include "grb.h"
// #include "DDSolver.h"

extern const Network network;

//

class OutObject{
public:
    double lb;
    double ub;
    vector<Node_t> nodes;
    bool success;

    OutObject(double lb_, double ub_, vector<Node_t> nodes_, bool s) :
        lb{lb_}, ub{ub_}, nodes{std::move(nodes_)}, success{s}{}
};

class NodeExplorer {

    //shared_ptr<Network> networkPtr;

    CutContainer feasibilityCuts;
    CutContainer optimalityCuts;
    GRBEnv env = GRBEnv();

    const shared_ptr<Network> networkPtr;

public:

    explicit NodeExplorer(const shared_ptr<Network>& networkPtr_) : networkPtr{networkPtr_}, feasibilityCuts{FEASIBILITY}, optimalityCuts{OPTIMALITY} {
        // env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag,0);
        env.set(GRB_IntParam_Threads,1);
    }

    NodeExplorer(const shared_ptr<Network>& networkPtr_, pair<CutContainer, CutContainer> cuts): networkPtr{networkPtr_}, feasibilityCuts{cuts.first}, optimalityCuts{cuts.second} {
        env.set(GRB_IntParam_OutputFlag,0);
        env.set(GRB_IntParam_Threads,1);
    }

	pair<CutContainer, CutContainer> getCuts() noexcept { return {feasibilityCuts, optimalityCuts};}


    OutObject process(Node_t node, double optimalLB);
    OutObject process2(Node_t node, double optimalLB);
    OutObject process3(Node_t node, double optimalLB);
    OutObject process4(Node_t node, double optimalLB);
    OutObject process5(Node_t node, double optimalLB);
    void clearCuts();

    #ifdef SOLVER_STATS
    void displayCutStats() const noexcept {
        cout << "Number of feasibility cuts in the container " << feasibilityCuts.cuts.size() << endl;
        cout << "Number of optimality cuts in the container: " << optimalityCuts.cuts.size() << endl;
    }
    #endif

};




#endif //NODEEXPLORER_H
