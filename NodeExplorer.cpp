//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include <complex>
#include "DD.h"


OutObject NodeExplorer::process(const Node_t node, double optimalLB) {

  	/*
		node should be eligible for processing.
  	 */

    double lowerBound = node.lb;
    double upperBound = node.ub;

    STEP_1: // refine relaxed tree with global feasibility cuts.
    {
        DDNode root1 {0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD1{networkPtr,EXACT};
        auto _ = relaxedDD1.build(root1,1);
        // refine tree with feasibility cuts
#ifdef DEBUG
        cout << "Applying Feasibility cuts heuristically on the relaxed tree." << endl;
#endif
        for (auto start = feasibilityCuts.cuts.rbegin(); start != feasibilityCuts.cuts.rend(); ++start) {
            const auto cut = *start;
        // for (const auto& cut: feasibilityCuts.cuts) {
            // cout << "-----------------------------------------------------------------------Zart" << endl;
            // if any of the cuts make the tree infeasible? get another node to explore.
            if (!relaxedDD1.applyFeasibilityCutHeuristic(cut)) {
#ifdef DEBUG
                cout << "Feasibility cut made tree infeasible" << endl;
#endif
                return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), {}, false}; // get next node
            }
        }
        // apply optimality cuts to the relaxed tree.
        DD relaxedDD2{networkPtr, EXACT};
        DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD2.build(root2,1);
        for(auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {

            upperBound = relaxedDD1.applyOptimalityCutHeuristic(*start);
            if (upperBound < optimalLB) {
                // break
#ifdef DEBUG
                cout << "upper bound < global lower bound." << endl;
#endif
                return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), {}, false};
            }
        }
    }

#ifdef DEBUG
    cout << "STEP_1 completed" << endl;
#endif
    STEP_2:
    DDNode root2 {0, node.globalLayer, node.states, node.solutionVector};
    // build restricted DD
    DD restrictedDD{networkPtr,RESTRICTED};
    auto cutset = restrictedDD.build(root2,1);
    // goto STEP_2C;
    // if(!cutset) goto STEP_2C;

    // cout << "STEP_2 completed" << endl;

    STEP_2A:
    // apply feasibility cuts on restricted DD.
    {
#ifdef DEBUG
        cout << "Applying feasibility cuts on the restricted tree in reverse order" << endl;
#endif
        // auto st = feasibilityCuts.cuts.rbegin();
        // auto end = (feasibilityCuts.cuts.size() > 100) ? feasibilityCuts.cuts.rbegin() + 100 : feasibilityCuts.cuts.rend();
        auto end = feasibilityCuts.cuts.rend();
        for (auto start = feasibilityCuts.cuts.rbegin(); start != end; ++start){
            const auto& cut = *start;
        // for (const auto& cut: feasibilityCuts.cuts) {
            if(!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                // tree is infeasible. a layer is removed. get cutset and set lb, ub to root's
                // auto cs = cutset.value();

                // for (auto& exactNode : cutset.value()) {
                //     exactNode.lb = node.lb;
                //     exactNode.ub = node.ub;
                // }
#ifdef DEBUG
                cout << "State < 0 after applying feasiblity cut on restricted DD" << endl;
#endif
                // return {node.lb, node.ub, cutset.value(), true};
                goto FINAL;
            }
        }
    }

    if (!cutset) goto STEP_2C; // if tree is exact, do not apply opt cuts on restricted tree.

    // cout << "STEP_2A completed" << endl;
    STEP_2B:
    {
#ifdef DEBUG
        cout << "Applying optimality cuts on restricted tree in reverse order." << endl;
#endif
        //applying optimality cuts in reverse order.
        // auto end = (optimalityCuts.cuts.size() > 100) ? optimalityCuts.cuts.rbegin()+100 : optimalityCuts.cuts.rend();
        auto end = optimalityCuts.cuts.rend();
        for (auto start = optimalityCuts.cuts.rbegin(); start != end; ++start) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start);
        }
    }
#ifdef DEBUG
    cout << "Lower bound after applying all optimality cuts: " << lowerBound<< endl;
#endif

    // // refinement
    STEP_2C:
    {
        vector<int> previousSolution;
#ifdef DEBUG
        cout << "Starting actual refinement." << endl;
#endif
        while (true) {
            // create relaxed tree and add refine with feasibility cuts
            //if (iter++ > 5) break;
            auto solution = restrictedDD.solution();
#ifdef DEBUG
            cout << "Solution: "; for (auto s: solution) cout << s << " "; cout << endl;
#endif
            if(previousSolution == solution) {
#ifdef DEBUG
                cout << "Previous solution returned" << endl;
#endif
                break;
            }
            previousSolution = solution;

            GuroSolver solver {networkPtr,env};
            auto y_bar = w2y(solution, networkPtr);
            auto cut = solver.solveSubProblemInstance(y_bar,0);
#ifdef DEBUG
            cout << "Cut: " << cut.cutType << " RHS: " << cut.RHS<< endl;
            for (const auto&[k,v] : cut.cutCoeff) {
                auto [i,q,j] = k;
                // cout << i << ", " << q << ", " << j << " : " << v << endl;
            }
#endif
            // auto cut = solver.solveSubProblem(y_bar);
            if (cut.cutType == FEASIBILITY) {
                // check if cut exists
                // if(feasibilityCuts.isCutExists(cut)) {
                //     break; // does it happen?
                // }
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
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
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                // LATER if lower bound is < global lower bound, break it
                optimalityCuts.insertCut(cut);
            }
        }
    }
#ifdef DEBUG
    cout << "Step 2 completed" << endl;
#endif
    STEP_3:
    {
        DDNode root3 {0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD2{networkPtr, EXACT};
        relaxedDD2.build(root3,1);
        //double upperBound = node.ub;
#ifdef DEBUG
        cout << "Applying optimality cuts on the relaxed tree heuristically." << endl;
#endif
        for (auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {
            upperBound = relaxedDD2.applyOptimalityCutHeuristic(*start);
        }
    }

#ifdef DEBUG
    cout << "lower bound : " << lowerBound << " , upper bound: " << upperBound << endl;
#endif
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


OutObject NodeExplorer::process2(const Node_t node, const double optimalLB) {

    // relaxed Tree 1,
    double lowerBound = node.lb;
    double upperBound = node.ub;

    STEP_1:

        DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD1{networkPtr, EXACT};
        relaxedDD1.build(root1,1);

        // apply cuts heuristically.
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();

        for (; start != end; ++start) {
            if (!relaxedDD1.applyFeasibilityCutHeuristic(*start)) {
#ifdef DEBUG
                cout << "tree became infeasible" << endl;
#endif
                return {0, 0, {}, false};
            }
        }


#ifdef DEBUG
    cout << "Step 1 completed" << endl;
#endif

    DD restrictedDD{networkPtr, RESTRICTED};
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    auto cutset = restrictedDD.build(root2,1);

    if (cutset) {
        // build relaxed DD 2
        DDNode root3{0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD2{networkPtr, EXACT};
        relaxedDD2.build(root3,1);

        // apply optimality cuts heuristically
        auto start = optimalityCuts.cuts.rbegin();
        auto end = optimalityCuts.cuts.rend();

        for (; start != end; ++start) {
            upperBound = relaxedDD2.applyOptimalityCutHeuristic(*start);
             if( upperBound < optimalLB) {
#ifdef DEBUG
                 cout << "upperbound : " << upperBound << " , < optimal LB: " << optimalLB << endl;
#endif
                 return {std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), {}, false};
             }
        }
    }

    STEP_2:
    {
        // apply cuts on restricted tree.

        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();

        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                lowerBound = node.lb;
                // upperBound = upperBound;
                goto FIN;
            }
        }

        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();

        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) goto FIN;
        }
    }
    // cout << "STEP 2 completed" << endl;

    STEP_3: // actual refinement
    {
        vi previousSolution;

        while (true) {

            vi solution = restrictedDD.solution();

            if (previousSolution == solution) {
#ifdef DEBUG
               cout << "previous solution returned" << endl;
#endif
                break;
            }
            previousSolution = solution;

            auto y_bar = w2y(solution, networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar, 0);

            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!relaxedDD1.applyFeasibilityCutHeuristic(cut)){return {0,0,{}, false};}

                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    // tree became infeasible.
                    lowerBound = node.lb;
                    // cout << "Feasibility cut made tree infeasible." << endl;
                    // upperBound = no;
                    goto FIN;
                }
            }
            else {

                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                optimalityCuts.insertCut(cut);
                if (lowerBound <= optimalLB) {
                    //
                    lowerBound = node.lb;
                    break;
                }
            }
        }
    }

    FIN:
    {
        // if cutset is not empty, update the lower boudn
        if (cutset) {
            for (auto& e : cutset.value()) {
                e.ub = upperBound;
                e.lb = lowerBound;
            }
            return {lowerBound, upperBound, cutset.value(), true};
        }
        return {lowerBound, upperBound, {}, true};
    }
    // apply cuts heuristically

}

OutObject NodeExplorer::process3(const Node_t node, const double optimalLB) {
    double lowerBound = node.lb;
    double upperBound = node.ub;

    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1,1);

    DD relaxedDD{networkPtr, EXACT};

    if (cutset) {
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2,1);

        // apply cuts on the relaxed DD. if any cut make DD infeasible, process new node.
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();

        for (; start != end; ++start) {
            if (!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
        }

        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();

        for (; start2 != end2; ++start2) {
            upperBound = relaxedDD.applyOptimalityCutHeuristic(*start2);
            // this node cannot do better.
            if (upperBound <= optimalLB) { return {0,0,{}, false};}
        }
    }

    // apply cuts to restricted tree.
    auto start = feasibilityCuts.cuts.rbegin();
    auto end = feasibilityCuts.cuts.rend();

    for (; start != end; ++start) {
        if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
            if (cutset) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    // c_node.ub = node.ub;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
            return{node.lb, upperBound, {}, true};
        }
    }

    auto start2 = optimalityCuts.cuts.rbegin();
    auto end2 = optimalityCuts.cuts.rend();

    for (; start2 != end2; ++start2) {
        lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
    }

    /// actual refinement ///

    vi previousSolution;

    while (true) {
        vi solution = restrictedDD.solution();
        if (previousSolution == solution) {break;}
        previousSolution = solution;

        auto y_bar = w2y(solution , networkPtr);
        GuroSolver solver{networkPtr, env};
        auto cut = solver.solveSubProblemInstance(y_bar, 0);

        if (cut.cutType == FEASIBILITY) {
            feasibilityCuts.insertCut(cut);
            if (cutset){if (!relaxedDD.applyFeasibilityCutHeuristic(cut)) return {0,0,{}, false};}
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                // copy lower bound from parent, upper bound achieved so far.
                if (cutset) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb, upperBound, cutset.value(), true};
                }
                return {node.lb,upperBound,{}, true};
            }

        }
        else {
            optimalityCuts.insertCut(cut);
            if (cutset) upperBound = relaxedDD.applyOptimalityCutHeuristic(cut);
            if (upperBound <= optimalLB) {
                return {0,0,{}, false};
            }
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
        }
    }

    if (cutset) {
        for (auto& c_node: cutset.value()) {
            c_node.ub = upperBound;
            c_node.lb = lowerBound;
        }
        return {lowerBound, upperBound, cutset.value(), true};
    }
    return {lowerBound, upperBound, {}, true};

}
void NodeExplorer::clearCuts() {
    feasibilityCuts.clearContainer();
    optimalityCuts.clearContainer();
}

