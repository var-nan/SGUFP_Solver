//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include <complex>
#include "DD.h"


OutObject NodeExplorer::process( const Network& network, const Node_t node) {

  	/*
		node should be eligible for processing.
  	 */

    double lowerBound = node.lb;
    double upperBound = node.ub;

    STEP_1: // refine relaxed tree with global feasibility cuts.
    DDNode root1 {0, node.globalLayer, node.states, node.solutionVector};
    DD relaxedDD1{EXACT};
    auto _ = relaxedDD1.build(network, root1);
    // refine tree with feasibility cuts
    for (auto start = feasibilityCuts.cuts.rbegin(); start != feasibilityCuts.cuts.rend(); ++start) {
        const auto cut = *start;
    // for (const auto& cut: feasibilityCuts.cuts) {
        // cout << "-----------------------------------------------------------------------Zart" << endl;
        // if any of the cuts make the tree infeasible? get another node to explore.
        if (!relaxedDD1.applyFeasibilityCutHeuristic(network, cut)) {
            return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), {}, false}; // get next node
        }
    }

    // cout << "STEP_1 completed" << endl;

    STEP_2:
    DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};
    // build restricted DD
    DD restrictedDD{RESTRICTED};
    auto cutset = restrictedDD.build(network, root2);
    // goto STEP_2C;
    // if(!cutset) goto STEP_2C;

    // cout << "STEP_2 completed" << endl;

    STEP_2A:
    // apply feasibility cuts on restricted DD.
    for (auto start = feasibilityCuts.cuts.rbegin(); start != feasibilityCuts.cuts.rend(); ++start){
        const auto cut = *start;
    // for (const auto& cut: feasibilityCuts.cuts) {
        if(!restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut)) {
            // tree is infeasible. a layer is removed. get cutset and set lb, ub to root's
            // auto cs = cutset.value();

            // for (auto& exactNode : cutset.value()) {
            //     exactNode.lb = node.lb;
            //     exactNode.ub = node.ub;
            // }
            // // cout << "State < 0 after applying feasiblity cut on restricted DD" << endl;
            // return {node.lb, node.ub, cutset.value(), true};
            goto FINAL;
        }
    }

    // cout << "STEP_2A completed" << endl;
    STEP_2B:
    {
        //applying optimality cuts in reverse order.
        for (auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(network, *start);
        }
    }
    // // for (const auto& cut: optimalityCuts.cuts) {
    // //     lowerBound= restrictedDD.applyOptimalityCutRestrictedLatest(network, cut);
    // // }
    // // cout << "STEP_2B completed" << endl;
    // // cout << "Lower Bound so far: " << lowerBound << endl;
    // // refinement
    STEP_2C:
    {
    vector<int> previousSolution;
    // cout << "Starting refinement" << endl;
    while (true) {
        // create relaxed tree and add refine with feasibility cuts
        //if (iter++ > 5) break;
        auto solution = restrictedDD.solution();
        if(previousSolution == solution) {
            break;
        }
        previousSolution = solution;
        // for (auto prev_sol : solutions) {
        //     if (prev_sol == solution) {
        //         // break from loop.
        //         cout<< "================= previous solution returned =====================" << endl;
        //         goto FINAL;
        //     }
        // }
        // solutions.push_back(solution);
        // cout << "Solution from tree: "; for (auto s : solution) cout << s << " " ; cout << endl;
        GuroSolver solver {env, static_cast<int>(network.n)};
        auto y_bar = w2y(solution, network);
        auto cut = solver.solveSubProblemInstance(network,y_bar, 0);

        if (cut.cutType == FEASIBILITY) {
            // check if cut exists
            // if(feasibilityCuts.isCutExists(cut)) {
            //     break; // does it happen?
            // }
            feasibilityCuts.insertCut(cut);
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(network, cut)) {
                // tree is removed, return the cutset and lower bound and upper bound.
                lowerBound  = node.lb;
                goto FINAL;
                // goto FINAL;
                // if (cutset) goto FINAL;
                //     //return {lowerBound, node.ub, cutset.value(), true};
                // return {lowerBound, node.lb, {}, true};
            }
        }
        else {
            // if (optimalityCuts.isCutExists(cut)) {
            //     // cout << "Optimality cut exists" << endl;
            //     break; // set lower bound
            // }
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(network, cut);
            // LATER if lower bound is < global lower bound, break it
            optimalityCuts.insertCut(cut);
        }
    }
    }
    STEP_3:
    {
        DDNode root3 {0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD2{EXACT};
        relaxedDD2.build(network, root3);
        //double upperBound = node.ub;
        for (auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {
            upperBound = relaxedDD2.applyOptimalityCutHeuristic(network, *start);
        }
    }
    // for (const auto& cut : optimalityCuts.cuts) {
    //     upperBound = relaxedDD2.applyOptimalityCutHeuristic(network, cut);
    //     // LATER if the upper bound is < global lower bound, then break it.
    // }
    // cout << "-----------------------------------Upper1 Bound so far: " << upperBound<< endl;


    FINAL:
    // optimalityCuts.clearContainer();
    //upperBound = node.ub;

    // cout << "One instance of node explorer finished" << endl;
    if (cutset) {
        // state is upper bound.
        for (auto& cutNode : cutset.value()) {
           cutNode.ub = upperBound;
            cutNode.lb = lowerBound;
        }
            return {lowerBound, upperBound, cutset.value(), true};
    }
    return  {lowerBound,upperBound, {}, true};
}

void NodeExplorer::clearCuts() {
    feasibilityCuts.clearContainer();
    optimalityCuts.clearContainer();
}

