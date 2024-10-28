//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include "DD.h"


void NodeExplorer::process( Node_t node) {

    DDNode root1 {0, node.globalLayer, node.states, node.solutionVector};
    CutContainer globalOptCuts{OPTIMALITY};
    CutContainer globalFeasCuts{FEASIBILITY};

    CutContainer optCuts{OPTIMALITY};
    CutContainer feasCuts{FEASIBILITY};

    double lowerBound;
    double upperBound;

    // phase 1: refine relaxed tree with global feasibility cuts.
    // TODO: ASAP create copy of root for multiple DD.
    DD relaxedDD1{RELAXED};
    relaxedDD1.build(network, root1);
    // refine tree
    bool inFeasible = false;
    for (const auto& cut: globalFeasCuts.cuts) {

        relaxedDD1.applyFeasibilityCutRestricted(network, cut);

        if (inFeasible) {
            // exit and get next node.
            return;
        }
    }

    DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};

    // build restricted DD
    DD restrictedDD{RESTRICTED};
    restrictedDD.build(network, root2);

    // refinement

    while (true) {
        // create relaxed tree and add refine with feasibility cuts

        auto solution = restrictedDD.solution();

        // const GRBEnv env1 = GRBEnv();
        // GRBModel mod = GRBModel(env1);


    }



    DDNode root3 {0, node.globalLayer, node.states, node.solutionVector};



}