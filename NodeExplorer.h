//
// Created by nandgate on 10/27/2024.
//

#ifndef NODEEXPLORER_H
#define NODEEXPLORER_H

#include "DD.h"
#include "Cut.h"
#include "gurobi_c++.h"
// #include "DDSolver.h"

extern const Network network;



class Pavani{
public:
    double lb;
    double ub;
    vector<Node_t> nodes;
    bool success;
};

class NodeExplorer {

    //shared_ptr<Network> networkPtr;

    CutContainer feasibilityCuts;
    CutContainer optimalityCuts;
    GRBEnv env = GRBEnv();

public:

    NodeExplorer() : feasibilityCuts{FEASIBILITY}, optimalityCuts{OPTIMALITY} {
        // env = GRBEnv();
    }

    Pavani process(const Network& network, Node_t node);

    void doSomething() {

    }

};




#endif //NODEEXPLORER_H
