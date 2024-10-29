//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include "DD.h"


Pavani NodeExplorer::process( const Network& network, Node_t node) {

  	/*
		node should be eligible for processing.
  	 */

    // check if restricted tree compiled with this node is exact?

    double lowerBound;
    double upperBound;

    // phase 1: refine relaxed tree with global feasibility cuts.
    // TODO: ASAP create copy of root for multiple DD.
    DDNode root1 {0, node.globalLayer, node.states, node.solutionVector};
    DD relaxedDD1{RELAXED};
    relaxedDD1.build(network, root1);
    // refine tree with feasibility cuts
    for (const auto& cut: feasibilityCuts.cuts) {
        if (relaxedDD1.applyFeasibilityCutHeuristic(network, cut)) return; // get next node
    }

    DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};
    // build restricted DD
    DD restrictedDD{RESTRICTED};
    restrictedDD.build(network, root2);

    // apply feasibility cuts on restricted DD.
    for (const auto& cut: feasibilityCuts.cuts) {
        restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut);
        // TODO handle case if the entire tree is removed or not.
    }

    // apply optimality cuts on restricted DD.
    for (const auto& cut: optimalityCuts.cuts) {
        restrictedDD.applyOptimalityCutRestrictedLatest(network, cut);
        // TODO handle case if the LB is less than the global LB, break and output cutset.
    }

    // refinement

    while (true) {
        // create relaxed tree and add refine with feasibility cuts

        auto solution = restrictedDD.solution();

        // const GRBEnv env1 = GRBEnv();
        // GRBModel mod = GRBModel(env1);


    }



    DDNode root3 {0, node.globalLayer, node.states, node.solutionVector};



}