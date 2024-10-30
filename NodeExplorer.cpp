//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include <complex>
#include "DD.h"


OutObject NodeExplorer::process( const Network& network, Node_t node) {

  	/*
		node should be eligible for processing.
  	 */

    // check if restricted tree compiled with this node is exact?

    cout << "Partial Solution of node to be processed: "; for (auto s : node.solutionVector) cout << s << " ";
    cout << endl;

    double lowerBound = node.lb;

    STEP_1: // refine relaxed tree with global feasibility cuts.
    DDNode root1 {0, node.globalLayer, node.states, node.solutionVector};
    DD relaxedDD1{EXACT};
    auto _ = relaxedDD1.build(network, root1);
    // refine tree with feasibility cuts
    for (const auto& cut: feasibilityCuts.cuts) {
        // if any of the cuts make the tree infeasible? get another node to explore.
        if (!relaxedDD1.applyFeasibilityCutHeuristic(network, cut)) {
            // cout << "State < 0 after applying feasibilty cut heuristically on relaxed DD. " << endl;
            return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), {}, false}; // get next node
        }
    }

    // cout << "STEP_1 completed" << endl;

    STEP_2:
    DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};
    // build restricted DD
    DD restrictedDD{RESTRICTED};
    auto cutset = restrictedDD.build(network, root2);

#ifdef DEBUG
    // if (cutset) {
    //     cout << "Node layer for this node " << cutset.value()[0].globalLayer << endl;
    //     cout << "partial solution of first node from cutset: ";
    //     for (auto s: cutset.value()[0].solutionVector) cout << s << " "; cout << endl;
    // }
#endif

    // cout << "STEP_2 completed" << endl;

    STEP_2A:
    // apply feasibility cuts on restricted DD.
    for (const auto& cut: feasibilityCuts.cuts) {
        if(!restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut)) {
            // tree is infeasible. a layer is removed. get cutset and set lb, ub to root's
            // auto cs = cutset.value();
            for (auto& exactNode : cutset.value()) {
                exactNode.lb = node.lb;
                exactNode.ub = node.ub;
            }
            // cout << "State < 0 after applying feasiblity cut on restricted DD" << endl;
            return {node.lb, node.ub, cutset.value(), false};
        }
    }

    // cout << "STEP_2A completed" << endl;
    STEP_2B:
    // apply optimality cuts on restricted DD.
    for (const auto& cut: optimalityCuts.cuts) {
        lowerBound= restrictedDD.applyOptimalityCutRestrictedLatest(network, cut);
    }
    // cout << "STEP_2B completed" << endl;
    // cout << "Lower Bound so far: " << lowerBound << endl;
    // refinement
    STEP_2C:
    int iter = 0;
    // cout << "Starting refinement" << endl;
    while (true) {
        // create relaxed tree and add refine with feasibility cuts
        if (iter++ > 5) break;
        auto solution = restrictedDD.solution();
        // cout << "Solution from tree: "; for (auto s : solution) cout << s << " " ; cout << endl;
        GuroSolver solver {env, static_cast<int>(network.n)};
        auto y_bar = w2y(solution, network);
        auto cut = solver.solveSubProblemInstance(network,y_bar, 0);
        cout << "RHS: " << cut.RHS << endl;
        if (cut.cutType == FEASIBILITY) cout << "Type: FEASIBILITY" << endl;
        else cout << "Type: OPTIMALITY" << endl;
        for (auto [k,v] : cut.cutCoeff) {
            auto [i,q,j] = k;
            cout << "i: " << i << ", q: " << q << " j: " << j << " val: " << v << endl;
        }
        if (cut.cutType == FEASIBILITY) {
            // check if cut exists
            if(feasibilityCuts.isCutExists(cut)) {
                break; // does it happen?
            }
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut)) {
                // tree is removed, return the cutset and lower bound and upper bound.
                return {lowerBound, node.ub, cutset.value(), true};
            }
            feasibilityCuts.insertCut(cut);
        }
        else {
            if (optimalityCuts.isCutExists(cut)) {
                // cout << "Optimality cut exists" << endl;
                break; // set lower bound
            }
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(network, cut);
            optimalityCuts.insertCut(cut);
        }
    }
    // cout << "Lower Bound so far: " << lowerBound << endl;
 // cout << "STEP_2C completed" << endl;
    STEP_3:
    DDNode root3 {0, node.globalLayer, node.states, node.solutionVector};
    DD relaxedDD2{EXACT};
    auto cutset3 = relaxedDD2.build(network, root3);
    double upperBound = node.ub;
    for (const auto& cut : optimalityCuts.cuts) {
        upperBound = relaxedDD2.applyOptimalityCutHeuristic(network, cut);
    }
    // cout << "-----------------------------------Upper1 Bound so far: " << upperBound<< endl;



    // cout << "One instance of node explorer finished" << endl;
    if (cutset) {
        // state is upper bound.
    for (auto& cutNode : cutset.value()) {
       cutNode.ub = upperBound;
    }
        return {lowerBound, upperBound, cutset.value(), true};
    }
    return  {lowerBound,upperBound, {}, false};
}

