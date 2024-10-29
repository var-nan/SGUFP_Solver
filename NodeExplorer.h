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

typedef struct Node {
    vi states;
    vi solutionVector;
    double lb;
    double ub;
    uint globalLayer;

    Node(){}

    Node(vi states_, vi solutionVector_, double lb_, double ub_, uint globalLayer_):
        states{std::move(states_)}, solutionVector{std::move(solutionVector_)},
        lb{lb_}, ub{ub_}, globalLayer{globalLayer_}{}
    // int a[40];
} Node_t;

class Pavani{
public:
    double lb;
    double ub;
    vector<Node_t> nodes;
    bool success;
};

class NodeExplorer {

    //shared_ptr<Network> networkPtr;

    CutContainer<Cut> feasibilityCuts;
    CutContainer<Cut> optimalityCuts;
    GRBEnv env;

public:

    NodeExplorer() : feasibilityCuts{FEASIBILITY}, optimalityCuts{OPTIMALITY} {
        env = GRBEnv();
    }

    Pavani process(const Network& network, Node_t node);

    void doSomething() {

    }

};




#endif //NODEEXPLORER_H
