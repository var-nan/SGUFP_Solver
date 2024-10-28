//
// Created by nandgate on 10/27/2024.
//

#ifndef NODEEXPLORER_H
#define NODEEXPLORER_H

#include "DD.h"
#include "Cut.h"
// #include "DDSolver.h"

extern const Network network;

typedef struct Node {
    vi states;
    vi solutionVector;
    double lb;
    double ub;
    uint globalLayer;
    // int a[40];
} Node_t;

class NodeExplorer {

    //shared_ptr<Network> networkPtr;

public:
    void process(Node_t node);

    void doSomething() {

    }

};




#endif //NODEEXPLORER_H
