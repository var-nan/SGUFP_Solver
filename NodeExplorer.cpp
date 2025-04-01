//
// Created by nandgate on 10/27/2024.
//

#include "NodeExplorer.h"
#include <complex>
#include "DD.h"
#include <algorithm>
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
        auto _ = relaxedDD1.build(root1);
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
        relaxedDD2.build(root2);
        for(auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {

            upperBound = relaxedDD1.applyOptimalityCutHeuristic(*start,optimalLB,upperBound);
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
    auto cutset = restrictedDD.build(root2);
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
            auto cut = solver.solveSubProblemInstance(y_bar);
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
        relaxedDD2.build(root3);
        //double upperBound = node.ub;
#ifdef DEBUG
        cout << "Applying optimality cuts on the relaxed tree heuristically." << endl;
#endif
        for (auto start = optimalityCuts.cuts.rbegin(); start != optimalityCuts.cuts.rend(); ++start) {
            upperBound = relaxedDD2.applyOptimalityCutHeuristic(*start,optimalLB,upperBound);
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
        relaxedDD1.build(root1);

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
    auto cutset = restrictedDD.build(root2);

    if (cutset) {
        // build relaxed DD 2
        DDNode root3{0, node.globalLayer, node.states, node.solutionVector};
        DD relaxedDD2{networkPtr, EXACT};
        relaxedDD2.build(root3);

        // apply optimality cuts heuristically
        auto start = optimalityCuts.cuts.rbegin();
        auto end = optimalityCuts.cuts.rend();

        for (; start != end; ++start) {
            upperBound = relaxedDD2.applyOptimalityCutHeuristic(*start ,optimalLB,upperBound);
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
            auto cut = solver.solveSubProblemInstance(y_bar);

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
    auto cutset = restrictedDD.build(root1);

    DD relaxedDD{networkPtr, EXACT};

    if (cutset) {
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);

        // apply cuts on the relaxed DD. if any cut make DD infeasible, process new node.
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();

        for (; start != end; ++start) {
            if (!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
        }

        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();

        for (; start2 != end2; ++start2) {
            upperBound = relaxedDD.applyOptimalityCutHeuristic(*start2 ,optimalLB,upperBound);
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
        auto cut = solver.solveSubProblemInstance(y_bar);

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
            if (cutset) upperBound = relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound);
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




OutObject NodeExplorer::process4(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ub: " << node.ub  << " cur opt: " << optimalLB << "        ";
    // cout << "process 4 start" << endl;
    double lowerBound = node.lb;
    double upperBound = node.ub;
    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1);
    if(restrictedDD.isTreeExact()) {
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) return {0,0,{}, false};
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) { return {node.lb,upperBound,{}, false};}
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    return {lowerBound,upperBound,{}, true};
                }
                previousLowerBound=lowerBound;
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    return {0,0,{}, true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) { return {node.lb,upperBound,{}, false};}
            }
        }
    }else {
        // cout << "went to else" << endl;
        DD relaxedDD{networkPtr, EXACT};
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);
        cout << relaxedDD.tree[relaxedDD.tree.size()-2].size() << "  " ;
        // applying old feas cuts
        // cout << endl;
        // for (auto item:relaxedDD.tree) {
        //     cout << item.size() << endl;
        // }
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
                }
                cout << " sec fes  ";
                auto start = feasibilityCuts.cuts.rbegin();
                auto end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
                }
                auto start2 = optimalityCuts.cuts.rbegin();
                auto end2 = optimalityCuts.cuts.rend();
                // cout << "here 2" << endl;
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
                }
                // cout << "here 3" << endl;
                auto fub = upperBound;
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                // cout << "here 2" << endl;
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
                }
                auto sub = upperBound;
                if (fub != sub) {
                    cout << " ****** " ;
                }
                // cout << " here 1 " << endl;
                vi solution = relaxedDD.solution();
                auto y_bar = w2y(solution , networkPtr);
                GuroSolver solver{networkPtr, env};
                auto cut = solver.solveSubProblemInstance(y_bar);
                // cout << " here 2 " << endl;
                if (cut.cutType == OPTIMALITY) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                    optimalityCuts.insertCut(cut);
                    cout << " got opt in relax ref ";
                    if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
                }else {
                    cout << " got  feas in relax ref ";
                    if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                        feasibilityCuts.insertCut(cut);
                        return {0,0,{}, false};
                    }
                }
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        cout << " sec fes  ";
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
        }

        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        // cout << "here 2" << endl;
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2 ,optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
        }
        auto fub = upperBound;
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        // cout << "here 2" << endl;
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2 ,optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
            // lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
        }
        auto sub = upperBound;
        if (fub != sub) {
            cout << " ****** " ;
        }
        // cout << "here 2" << endl;
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {lowerBound,upperBound,cutset.value(), true};
                }else{
                    previousLowerBound=lowerBound;
                }
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {0,0,{}, false};
                }
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb,upperBound,cutset.value(), true};
                    // return {node.lb,upperBound,cutset.value(), true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB ) {return {node.lb,upperBound,{}, false};}
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
            }
        }
    }
}

OutObject NodeExplorer::process5(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";

    double lowerBound = node.lb;
    double upperBound = node.ub;

    // GuroSolver solver{networkPtr, env};
    // upperBound = solver.solvepartialProblem(node.solutionVector ,*networkPtr);
    // if (upperBound <= optimalLB) {
    //     return {node.lb,upperBound,{}, false};
    // }

    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1);

    // for (auto item : restrictedDD.tree) {
    //     cout << item.size() << endl;
    // }
    if(restrictedDD.isTreeExact()) {
        // cout << "exact tree detected!" << endl;
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                // cout << "exact tree pruned by feas!" << endl;
                return {0,0,{}, false};
            }
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) {
                // cout << "exact tree pruned by opt!" << endl;
                return {0,0,{}, false};
            }
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    // cout << "exact tree found a new solution" << endl;
                    return {lowerBound,upperBound,{}, true};
                }
                previousLowerBound=lowerBound;
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    // cout << "exact tree pruned by feas actual ref!" << endl;
                    return {0,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) {
                    // cout << "exact tree pruned by opt actual ref!" << endl;
                    return {0,0,{}, false};
                }
            }
        }
    }else {
        // cout << "went to else" << endl;
        DD relaxedDD{networkPtr, EXACT};
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);

        // cout<<endl;
        // cout << "relaxed DD"<< endl;
        // for (auto item: relaxedDD.tree) {
        //     cout << item.size() << endl;
        // }
        // cout << endl;

        cout << " " << relaxedDD.tree[relaxedDD.tree.size()-2].size() << " " ;
        /////////////////////////////////////////////////////////////////////////////////////
        /// mixing cut ///
        /// feas Cuts
        CutCoefficients Y_bar_coef;
        double rhs = 0.0;
        CutType type = FEASIBILITY;
        Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
        auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
        if (feasibilityCuts.cuts.size()>1) {
            specialFeasCut  = feasibilityCuts.cuts[0];
            specialFeasCut.RHS = 0.0;
            for(auto item : specialFeasCut.cutCoeff) {
                specialFeasCut.cutCoeff[item.first] = 0.0;
            }
            for (auto oneCut : feasibilityCuts.cuts) {
                specialFeasCut.RHS += oneCut.RHS* 1.0 ;
                for(auto item : oneCut.cutCoeff) {
                    specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
                }
            }
            // cout << " app special feas: " ;
            if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
                return {0,0,{}, false};
            }
        }

        /// opt Cuts
        CutCoefficients Y_bar_coef1;
        double rhs1 = 0.0;
        CutType type1 = OPTIMALITY ;
        Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
        auto TotaloptCutsNum = optimalityCuts.cuts.size();
        if (optimalityCuts.cuts.size()>1) {
            specialOptCut  = optimalityCuts.cuts[0];
            specialOptCut.RHS = 0.0;
            for(auto item : specialOptCut.cutCoeff) {
                specialOptCut.cutCoeff[item.first] = 0.0;
            }
            for (auto oneCut : optimalityCuts.cuts) {
                specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
                for(auto item : oneCut.cutCoeff) {
                    specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
                }
            }
            // cout << " app special opt: " ;
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        }
        /// mix of special opt and special feas
        CutCoefficients Y_bar_coef2;
        double rhs2 = 0;
        CutType type2 = OPTIMALITY;
        Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
        if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
            verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
            for(auto item : specialFeasCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] = 0.0;
            }
            for(auto item : specialFeasCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
            }
            for(auto item : specialOptCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
            }
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        }


        /////////////////////////////////////////////////////////////////////////////////////
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {0,0,{}, false};
            }
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {0,0,{}, false};
                    }
                }
                for (auto ct1 = 0; ct1 < 1 ; ++ct1) {
                    auto start = feasibilityCuts.cuts.rbegin();
                    auto end = feasibilityCuts.cuts.rend();
                    cout << " sfeas ";
                    for (; start != end; ++start) {
                        if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                            return {0,0,{}, false};
                        }
                    }
                }
                auto start2 = optimalityCuts.cuts.rbegin();
                auto end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                }
                auto start = feasibilityCuts.cuts.rbegin();
                auto end = feasibilityCuts.cuts.rend();
                // cout << " sfeas ";
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {0,0,{}, false};
                    }
                }
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                }
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {0,0,{}, false};
            }
        }
        for (auto ct1 = 0; ct1 < 1 ; ++ct1) {
            start = feasibilityCuts.cuts.rbegin();
            end = feasibilityCuts.cuts.rend();
            cout << " sfeas ";
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {0,0,{}, false};
                }
            }
        }
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        // cout << " sfeas ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {0,0,{}, false};
            }
        }
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        // cout << " sfeas ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {0,0,{}, false};
            }
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {lowerBound,upperBound,cutset.value(), true};
                }else{
                    previousLowerBound=lowerBound;
                }
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    // cout << " relaxed tree became infeas " ;
                    return {0,0,{}, false};
                }
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    // cout << " relaxed tree became infeas " ;
                    return {0,0,{}, false};
                }
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb,upperBound,cutset.value(), true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound),upperBound);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {
                    // cout << " relaxed pruned by opt! ";
                    return {0,0,{}, false};
                }
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    // cout << " restricted tree went below opt! " ;
                    return {node.lb, upperBound, cutset.value(), true};
                }
            }
        }
    }
}



OutObject NodeExplorer::process6(const Node_t node, const double optimalLB) {
    double lowerBound = node.lb;
    double upperBound = node.ub;
    DD relaxedDD{networkPtr, EXACT};
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    relaxedDD.build(root2);

    for (auto item : relaxedDD.tree) {
        cout << item.size() << endl;
    }
    // actual refinement
    double previousLowerBound;
    vi previousSolution;
    while (true) {
        vi solution = relaxedDD.solution();
        if (previousSolution == solution) {
            if (previousLowerBound == lowerBound) {
                return {lowerBound,upperBound,{}, true};
            }else{
                previousLowerBound=lowerBound;
            }
        }else {
            previousSolution = solution;
            previousLowerBound = lowerBound;
        }
        auto y_bar = w2y(solution , networkPtr);
        GuroSolver solver{networkPtr, env};
        auto cut = solver.solveSubProblemInstance(y_bar);
        if (cut.cutType == FEASIBILITY) {
            feasibilityCuts.insertCut(cut);
            if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                return {0,0,{}, false};
            }
            if (!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                return {node.lb,upperBound,{}, true};
            }
        }else {
            optimalityCuts.insertCut(cut);
            upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);;
            if (upperBound <= optimalLB) {return {0,0,{}, false};}
            lowerBound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
            cout << lowerBound << endl;
        }
    }
}
// from 4
OutObject NodeExplorer::process7(const Node_t node, const double optimalLB) {
    cout << "process 7 start" << endl;
    double lowerBound = node.lb;
    double upperBound = node.ub;
    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1);
    // cout << "passed restricted tree" << endl;
    if(restrictedDD.isTreeExact()) {
        // cout << "this is an exact tree" << endl;
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) return {0,0,{}, false};
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) { return {0,0,{}, false};}
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    return {lowerBound,upperBound,{}, true};
                }
                previousLowerBound=lowerBound;
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    return {-100000000,100000000,{}, true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) { return {0,0,{}, false};}
            }
        }
    }else {
        // cout << "went to else" << endl;
        DD relaxedDD{networkPtr, EXACT};
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {
                cout << "pruned" << endl;
                return {0,0,{}, false};
            }
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    for (int erc = 0 ; erc < 10 ; erc++) {
                        solution = relaxedDD.solution();
                        auto y_bar = w2y(solution , networkPtr);
                        GuroSolver solver{networkPtr, env};
                        auto cut = solver.solveSubProblemInstance(y_bar);
                        if (cut.cutType == OPTIMALITY) {
                            auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                            if (newuppperbound < upperBound) {
                                upperBound = newuppperbound;
                                cout << "new upper bound idea reduced the previous upperbound!" << endl;
                                if (upperBound <= optimalLB) {
                                    cout << "new idea pruned a node " << endl;
                                    return {0,0,{}, false};
                                }
                            }
                        }else{break;}
                    }
                    // extra relaxed refinement
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = lowerBound;
                    }
                    return {lowerBound,upperBound,cutset.value(), true};
                }else{
                    previousLowerBound=lowerBound;
                }
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {0,0,{}, false};
                }
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb,upperBound,cutset.value(), true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {0,0,{}, false};}
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
            }
        }
    }
}

OutObject NodeExplorer::process8(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ub: " << node.ub  << "    ";
    // cout << "process 8 started" << endl;
    // cout << "Layer: " << node.globalLayer << endl;
    double lowerBound = node.lb;
    double upperBound = node.ub;
    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1);
    if(restrictedDD.isTreeExact()) {
        cout << "this is an exact tree" << endl;
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) return {0,0,{}, false};
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) { return {0,0,{}, false};}
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    return {lowerBound,upperBound,{}, true};
                }
                previousLowerBound=lowerBound;
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    return {-100000000,100000000,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) { return {0,0,{}, false};}
            }
        }
    }else {
        // cout << "went to else" << endl;
        DD relaxedDD{networkPtr, EXACT};
        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);
        // for(auto item : relaxedDD.tree) {
        //     cout << item.size() << endl;
        // }
        // cout << endl;
        // cout << "end of building relaxed DD" << endl;
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) return {0,0,{}, false};
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                auto start2 = optimalityCuts.cuts.rbegin();
                auto end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                }

                for (int erc = 0 ; erc < node.globalLayer ; erc++) {
                    auto solution = relaxedDD.solution();
                    auto y_bar = w2y(solution , networkPtr);
                    GuroSolver solver{networkPtr, env};
                    auto cut = solver.solveSubProblemInstance(y_bar);
                    if (cut.cutType == OPTIMALITY) {
                        auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                        if (newuppperbound < upperBound) {
                            optimalityCuts.insertCut(cut);
                            cout << "new upper bound idea reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound << endl;
                            upperBound = newuppperbound;
                            if (upperBound <= optimalLB) {
                                cout << "new idea pruned a node " << endl;
                                return {0,0,{}, false};
                            }
                        }
                    }else {
                        // cout << "got Feasssssssss! " << endl;
                        break;
                    }
                }


                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }

        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {0,0,{}, false};}
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);

            if (lowerBound <= optimalLB) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    // extra relaxed refinement
                    for (int erc = 0 ; erc < node.globalLayer ; erc++) {
                        solution = relaxedDD.solution();
                        auto y_bar = w2y(solution , networkPtr);
                        GuroSolver solver{networkPtr, env};
                        auto cut = solver.solveSubProblemInstance(y_bar);
                        if (cut.cutType == OPTIMALITY) {
                            auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                            if (newuppperbound < upperBound) {
                                optimalityCuts.insertCut(cut);
                                cout << "new upper bound idea reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound << endl;
                                upperBound = newuppperbound;
                                if (upperBound <= optimalLB) {
                                    cout << "new idea pruned a node " << endl;
                                    return {0,0,{}, false};
                                }
                            }
                        }else {
                            // cout << "got Feasssssssss! " << endl;
                            break;
                        }
                    }
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = lowerBound;
                    }
                    return {lowerBound,upperBound,cutset.value(), true};
                }else{
                    previousLowerBound=lowerBound;
                }
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    // cout << "node got pruned by feas" << endl;
                    return {0,0,{}, false};
                }
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb,upperBound,cutset.value(), true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {
                    // cout << "node got pruned by ub"<< endl;
                    return {0,0,{}, false};
                }
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) {
                    // extra relaxed refinement
                    for (int erc = 0 ; erc < node.globalLayer ; erc++) {
                        solution = relaxedDD.solution();
                        auto y_bar = w2y(solution , networkPtr);
                        GuroSolver solver{networkPtr, env};
                        auto cut = solver.solveSubProblemInstance(y_bar);
                        if (cut.cutType == OPTIMALITY) {
                            auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound);
                            if (newuppperbound < upperBound) {
                                optimalityCuts.insertCut(cut);
                                cout << "new upper bound idea reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound << endl;
                                upperBound = newuppperbound;
                                if (upperBound <= optimalLB) {
                                    cout << "new idea pruned a node " << endl;
                                    return {0,0,{}, false};
                                }
                            }
                        }else {
                            // cout << " got Feassssssssss!   " << endl;
                            break;
                        }
                    }
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb, upperBound, cutset.value(), true};
                }
            }
        }
    }
}
OutObject NodeExplorer::process10(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer <<" ub: " << node.ub  << "             ";
    double lowerBound = node.lb;
    double upperBound = node.ub;
    DDNode root1{0, node.globalLayer, node.states, node.solutionVector};
    DD restrictedDD {networkPtr, RESTRICTED};
    auto cutset = restrictedDD.build(root1);
    // for (auto item : restrictedDD.tree) {
    //     cout << item.size() << endl;
    // }
    if(restrictedDD.isTreeExact()) {
        // cout << "exact tree detected!" << endl;
        // applying old feas cuts
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                // cout << "exact tree pruned by feas!" << endl;
                return {0,0,{}, false};
            }
        }
        // applying old opt cuts
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) {
                // cout << "exact tree pruned by opt!" << endl;
                return {0,0,{}, false};
            }
        }
        // actual refinement
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    // cout << "exact tree found a new solution" << endl;
                    return {lowerBound,upperBound,{}, true};
                }
                previousLowerBound=lowerBound;
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    // cout << "exact tree pruned by feas actual ref!" << endl;
                    return {0,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) {
                    // cout << "exact tree pruned by opt actual ref!" << endl;
                    return {0,0,{}, false};
                }
            }
        }
    }else {
        // cout << "went to else" << endl;
        DD relaxedDD{networkPtr, EXACT};

        DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
        relaxedDD.build(root2);

        cout << " " << relaxedDD.tree[relaxedDD.tree.size()-2].size() << " " ;
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            // cout << "here 1 " << endl;
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {0,0,{}, false};
            }
            // cout << "here 2 " << endl;
            if (!restrictedDD.applyFeasibilityCutRestrictedLatest(*start)) {
                auto start2 = optimalityCuts.cuts.rbegin();
                auto end2 = optimalityCuts.cuts.rend();
                auto prvUB = upperBound;
                auto lastEffectiveCut = start2;
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if(prvUB != upperBound) {
                        prvUB = upperBound;
                        lastEffectiveCut = start2;
                    }
                    if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                }
                if (optimalityCuts.cuts.size()>0 && node.globalLayer > 20) {
                    relaxedDD.applyOptimalityCutHeuristic(*lastEffectiveCut, optimalLB,upperBound);
                    for (int erc = 0 ; erc < 10 ; erc++) {
                        auto solution = relaxedDD.solution();
                        auto y_bar = w2y(solution , networkPtr);
                        GuroSolver solver{networkPtr, env};
                        auto cut = solver.solveSubProblemInstance(y_bar);
                        if (cut.cutType == OPTIMALITY) {
                            auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                            if (newuppperbound < upperBound) {
                                optimalityCuts.insertCut(cut);
                                cout << "reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound ;
                                upperBound = newuppperbound;
                                if (upperBound <= optimalLB) {
                                    cout << "new idea pruned a node " ;
                                    return {0,0,{}, false};
                                }
                            }
                        }else {
                            if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                                feasibilityCuts.insertCut(cut);
                                cout << "new idea pruned a node " ;
                                return {0,0,{}, false};
                            }
                        }
                    }
                }
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                // cout << " restricted tree became infeas! " ;
                return {node.lb, upperBound, cutset.value(), true};
            }
        }

        // applying old opt cuts
        // cout << "applying old cuts after feas" << endl;
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        auto prvUB = upperBound;
        auto lastEffectiveCut = start2;
        // cout << "zart 1" << endl;
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if(prvUB != upperBound) {
                prvUB = upperBound;
                lastEffectiveCut = start2;
            }
        }
        // cout << "zart 2" << endl;
        if (optimalityCuts.cuts.size()>0 && node.globalLayer > 20) {
            relaxedDD.applyOptimalityCutHeuristic(*lastEffectiveCut, optimalLB,upperBound);
            for (int erc = 0 ; erc < 10 ; erc++) {
                auto solution = relaxedDD.solution();
                auto y_bar = w2y(solution , networkPtr);
                GuroSolver solver{networkPtr, env};
                auto cut = solver.solveSubProblemInstance(y_bar);
                if (cut.cutType == OPTIMALITY) {
                    auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                    if (newuppperbound < upperBound) {
                        optimalityCuts.insertCut(cut);
                        cout << "reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound ;
                        upperBound = newuppperbound;
                        if (upperBound <= optimalLB) {
                            cout << "new idea pruned a node ";
                            return {0,0,{}, false};
                        }
                    }
                }else {
                    // cout << "here 4 "<< endl;
                    if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                        feasibilityCuts.insertCut(cut);
                        cout << "new idea pruned a node " ;
                        return {0,0,{}, false};
                    }
                }
            }
        }
        for (auto& c_node: cutset.value()) {
            c_node.ub = upperBound;
            c_node.lb = node.lb;
        }
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(*start2);
            if (lowerBound <= optimalLB) {
                for (auto& c_node: cutset.value()) {
                    c_node.ub = upperBound;
                    c_node.lb = node.lb;
                }
                return {node.lb, upperBound, cutset.value(), true};
            }
        }
        // actual refinement
        // cout << "actual refinement" << endl;
        double previousLowerBound;
        vi previousSolution;
        while (true) {
            vi solution = restrictedDD.solution();
            if (previousSolution == solution) {
                if (previousLowerBound == lowerBound) {
                    for (int erc = 0 ; erc < 10 ; erc++) {
                        solution = relaxedDD.solution();
                        auto y_bar = w2y(solution , networkPtr);
                        GuroSolver solver{networkPtr, env};
                        auto cut = solver.solveSubProblemInstance(y_bar);
                        if (cut.cutType == OPTIMALITY) {
                            auto newuppperbound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                            if (newuppperbound < upperBound) {
                                optimalityCuts.insertCut(cut);
                                cout << "reduced the previous upperbound in iter: " << erc << "  results  " << upperBound << " --> " << newuppperbound ;
                                upperBound = newuppperbound;
                                if (upperBound <= optimalLB) {
                                    cout << "new idea pruned a node ";
                                    return {0,0,{}, false};
                                }
                            }
                        }else {
                            if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                                feasibilityCuts.insertCut(cut);
                                cout << "new idea pruned a node ";
                                return {0,0,{}, false};
                            }
                        }
                    }
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {lowerBound,upperBound,cutset.value(), true};
                }else{
                    previousLowerBound=lowerBound;
                }
            }else {
                previousSolution = solution;
                previousLowerBound = lowerBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    // cout << " relaxed tree became infeas " ;
                    return {0,0,{}, false};
                }
                if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    return {node.lb,upperBound,cutset.value(), true};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut, optimalLB,upperBound),upperBound);;
                if (upperBound <= optimalLB) {
                    // cout << " relaxed pruned by opt! ";
                    return {0,0,{}, false};
                }
                lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
                if (lowerBound <= optimalLB) {
                    for (auto& c_node: cutset.value()) {
                        c_node.ub = upperBound;
                        c_node.lb = node.lb;
                    }
                    // cout << " restricted tree went below opt! " ;
                    return {node.lb, upperBound, cutset.value(), true};
                }
            }
        }
    }
}




OutObject NodeExplorer::processX1(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";
    // node.lb = std::numeric_limits<int64_t>::min();
    double lowerBound = std::numeric_limits<int64_t>::min();
    double upperBound = std::numeric_limits<int64_t>::max();
    DD relaxedDD{networkPtr, EXACT};
    // for(auto item : node.states) {
    //    cout << item << "~";
    // }

    // cout << "node.globalLayer" << node.globalLayer << endl;
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    relaxedDD.build(root2);

    // cout << " " << relaxedDD.tree[relaxedDD.tree.size()-2].size() << " " ;
    // for (auto item : relaxedDD.tree) {
    //     cout << item.size() << endl;
    // }
    // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
    //     cout << i << " - " << relaxedDD.tree[i].size() << endl;
    // }
    CutCoefficients Y_bar_coef3;
    double rhs3 = 0.0;
    CutType type3 = OPTIMALITY;
    Cut lastEffectiveOptCut = Cut{type3, rhs3, Y_bar_coef3};


    if (relaxedDD.RelaxedisExact) {
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {lowerBound,0,{}, false};
            }
        }
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        cout << " sfeas ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {lowerBound,0,{}, false};
            }
        }
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }
        /// just in case the given opt sol is not actually opt
        double previousUpperBound;
        vi previousSolution;
        while (true) {
            vi solution = relaxedDD.solution();
            if (previousSolution == solution) {
                if (previousUpperBound == upperBound) {
                    return {upperBound,upperBound,{}, true};
                }else{
                    previousUpperBound=upperBound;
                }
            }else {
                previousSolution = solution;
                previousUpperBound = upperBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
            }
        }
    }else {
        // cout << " passed1 " ;
        cout << " " << relaxedDD.tree[relaxedDD.tree.size()-2].size() << " " ;

        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        // cout << "before fs1 ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        // cout << "after fs1";
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        cout << " sf ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            auto newUB = relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound);
            if (newUB < upperBound ) {
                upperBound = newUB;
                lastEffectiveOptCut = *start2;
            }
            if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
        }
        cout << " so ";
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            auto newUB = relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound);
            if (newUB < upperBound ) {
                upperBound = newUB;
                lastEffectiveOptCut = *start2;
            }
            if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
        }

    }

    double previousUpperBound;
    vi previousSolution;
    vector<vi> all_solutions={};




    // cout << " passed2 " ;
    while (true) {
        // auto start2 = optimalityCuts.cuts.rbegin();
        // auto end2 = optimalityCuts.cuts.rend();
        // for (; start2 != end2; ++start2) {
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        // }
        vi solution = relaxedDD.solution();
        for (auto item : all_solutions) {
            if (item == solution && previousUpperBound == upperBound) {
                auto start = feasibilityCuts.cuts.rbegin();
                auto end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {std::numeric_limits<int64_t>::min(),0,{}, false};
                    }
                }
                auto start2 = optimalityCuts.cuts.rbegin();
                auto end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                }
                start = feasibilityCuts.cuts.rbegin();
                end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {std::numeric_limits<int64_t>::min(),0,{}, false};
                    }
                }
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                }
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                // CutCoefficients Y_bar_coef;
                // double rhs = 0.0;
                // CutType type = FEASIBILITY;
                // Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
                // auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
                // if (feasibilityCuts.cuts.size()>1) {
                //     specialFeasCut  = feasibilityCuts.cuts[0];
                //     specialFeasCut.RHS = 0.0;
                //     for(auto item : specialFeasCut.cutCoeff) {
                //         specialFeasCut.cutCoeff[item.first] = 0.0;
                //     }
                //     for (auto oneCut : feasibilityCuts.cuts) {
                //         specialFeasCut.RHS += oneCut.RHS* 1.0 ;
                //         for(auto item : oneCut.cutCoeff) {
                //             specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
                //         }
                //     }
                //     // cout << " app special feas: " ;
                //     if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
                //         return {0,0,{}, false};
                //     }
                // }
                //
                // /// opt Cuts
                // CutCoefficients Y_bar_coef1;
                // double rhs1 = 0.0;
                // CutType type1 = OPTIMALITY ;
                // Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
                // auto TotaloptCutsNum = optimalityCuts.cuts.size();
                // if (optimalityCuts.cuts.size()>1) {
                //     specialOptCut  = optimalityCuts.cuts[0];
                //     specialOptCut.RHS = 0.0;
                //     for(auto item : specialOptCut.cutCoeff) {
                //         specialOptCut.cutCoeff[item.first] = 0.0;
                //     }
                //     for (auto oneCut : optimalityCuts.cuts) {
                //         specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
                //         for(auto item : oneCut.cutCoeff) {
                //             specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
                //         }
                //     }
                //     // cout << " app special opt: " ;
                //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
                //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                // }
                // /// mix of special opt and special feas
                // CutCoefficients Y_bar_coef2;
                // double rhs2 = 0;
                // CutType type2 = OPTIMALITY;
                // Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
                // if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
                //     verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
                //     for(auto item : specialFeasCut.cutCoeff) {
                //         verySpecialCut.cutCoeff[item.first] = 0.0;
                //     }
                //     for(auto item : specialFeasCut.cutCoeff) {
                //         verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
                //     }
                //     for(auto item : specialOptCut.cutCoeff) {
                //         verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
                //     }
                //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
                //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
                // }
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                /////////////////////////////////////////////////////
                start = feasibilityCuts.cuts.rbegin();
                end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {std::numeric_limits<int64_t>::min(),0,{}, false};
                    }
                }
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                }
                start = feasibilityCuts.cuts.rbegin();
                end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {std::numeric_limits<int64_t>::min(),0,{}, false};
                    }
                }
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                }

                auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
                for (auto& c_node: cutset) {
                    c_node.ub = upperBound;
                    c_node.lb = optimalLB;
                }
                return {lowerBound,upperBound,cutset, true};
            }
        }
        all_solutions.push_back(solution);
        previousUpperBound = upperBound;
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        // CutCoefficients Y_bar_coef;
        // double rhs = 0.0;
        // CutType type = FEASIBILITY;
        // Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
        // auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
        // if (feasibilityCuts.cuts.size()>1) {
        //     specialFeasCut  = feasibilityCuts.cuts[0];
        //     specialFeasCut.RHS = 0.0;
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         specialFeasCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for (auto oneCut : feasibilityCuts.cuts) {
        //         specialFeasCut.RHS += oneCut.RHS* 1.0 ;
        //         for(auto item : oneCut.cutCoeff) {
        //             specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
        //         }
        //     }
        //     // cout << " app special feas: " ;
        //     if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
        //         return {0,0,{}, false};
        //     }
        // }
        //
        // /// opt Cuts
        // CutCoefficients Y_bar_coef1;
        // double rhs1 = 0.0;
        // CutType type1 = OPTIMALITY ;
        // Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
        // auto TotaloptCutsNum = optimalityCuts.cuts.size();
        // if (optimalityCuts.cuts.size()>1) {
        //     specialOptCut  = optimalityCuts.cuts[0];
        //     specialOptCut.RHS = 0.0;
        //     for(auto item : specialOptCut.cutCoeff) {
        //         specialOptCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for (auto oneCut : optimalityCuts.cuts) {
        //         specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
        //         for(auto item : oneCut.cutCoeff) {
        //             specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
        //         }
        //     }
        //     // cout << " app special opt: " ;
        //     // upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
        //     // if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        // }
        // /// mix of special opt and special feas
        // CutCoefficients Y_bar_coef2;
        // double rhs2 = 0;
        // CutType type2 = OPTIMALITY;
        // Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
        // if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
        //     verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
        //     }
        //     for(auto item : specialOptCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
        //     }
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        // }
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////
        /////////////////////////////////////////////////////


        auto y_bar = w2y(solution , networkPtr);
        GuroSolver solver{networkPtr, env};
        auto cut = solver.solveSubProblemInstance(y_bar);
        if (cut.cutType == FEASIBILITY) {
            feasibilityCuts.insertCut(cut);
            if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                return {lowerBound,0,{}, false};
            }
            if (!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                return {lowerBound,0,{}, true};
            }
            if (optimalityCuts.cuts.size()>0) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(lastEffectiveOptCut, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
            }
        }else {
            optimalityCuts.insertCut(cut);
            upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);;
            lastEffectiveOptCut = cut;
            if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
        }
    }
}



void NodeExplorer::clearCuts() {
    feasibilityCuts.clearContainer();
    optimalityCuts.clearContainer();
}

OutObject NodeExplorer::processX2(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";
    // node.lb = std::numeric_limits<int64_t>::min();
    double lowerBound = std::numeric_limits<int64_t>::min();
    double upperBound = std::numeric_limits<int64_t>::max();
    DD relaxedDD{networkPtr, EXACT};


    // cout << "node.globalLayer" << node.globalLayer << endl;
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    relaxedDD.build(root2);


    // for (auto item : relaxedDD.tree) {
    //     cout << item.size() << endl;
    // }
    // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
    //     cout << i << " - " << relaxedDD.tree[i].size() << endl;
    // }

    CutCoefficients Y_bar_coef3;
    double rhs3 = 0.0;
    CutType type3 = OPTIMALITY;
    Cut lastEffectiveOptCut = Cut{type3, rhs3, Y_bar_coef3};


    if (relaxedDD.RelaxedisExact) {
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {lowerBound,0,{}, false};
            }
        }
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        cout << " sfeas ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyOptimalityCutRestrictedLatest(*start)) {
                return {lowerBound,0,{}, false};
            }
        }
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }
        /// just in case the given opt sol is not actually opt
        double previousUpperBound;
        vi previousSolution;
        while (true) {
            vi solution = relaxedDD.solution();
            if (previousSolution == solution) {
                if (previousUpperBound == upperBound) {
                    return {upperBound,upperBound,{}, true};
                }else{
                    previousUpperBound=upperBound;
                }
            }else {
                previousSolution = solution;
                previousUpperBound = upperBound;
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
            }
        }
    }else {
        // cout << " passed1 " ;
        cout << " " << relaxedDD.tree[relaxedDD.tree.size()-2].size() << " " ;

        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        // cout << "before fs1 ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        // cout << "after fs1";
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        cout << " sf ";
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            auto newUB = relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound);
            if (newUB < upperBound ) {
                upperBound = newUB;
                lastEffectiveOptCut = *start2;
            }
            if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
        }
        cout << " so ";
        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            auto newUB = relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound);
            if (newUB < upperBound ) {
                upperBound = newUB;
                lastEffectiveOptCut = *start2;
            }
            if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
        }

    }

    double previousUpperBound;
    vi previousSolution;
    vector<vi> all_solutions={};




    // cout << " passed2 " ;
    while (true) {
        // auto start2 = optimalityCuts.cuts.rbegin();
        // auto end2 = optimalityCuts.cuts.rend();
        // for (; start2 != end2; ++start2) {
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        // }

        // auto start = feasibilityCuts.cuts.rbegin();
        // auto end = feasibilityCuts.cuts.rend();
        // for (; start != end; ++start) {
        //     if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
        //         return {std::numeric_limits<int64_t>::min(),0,{}, false};
        //     }
        // }

        // CutCoefficients Y_bar_coef1;
        // double rhs1 = 0.0;
        // CutType type1 = OPTIMALITY ;
        // Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
        // auto TotaloptCutsNum = optimalityCuts.cuts.size();
        // if (optimalityCuts.cuts.size()>1) {
        //     specialOptCut  = optimalityCuts.cuts[0];
        //     specialOptCut.RHS = 0.0;
        //     for(auto item : specialOptCut.cutCoeff) {
        //         specialOptCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for (auto oneCut : optimalityCuts.cuts) {
        //         specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
        //         for(auto item : oneCut.cutCoeff) {
        //             specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
        //         }
        //     }
        //     // cout << " app special opt: " ;
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        // }
        // relaxedDD.applyOptimalityCutHeuristic(lastEffectiveOptCut, optimalLB,upperBound);
        auto solutions = relaxedDD.solution2();
        vector<vi> uniqueSols= {};
        for (auto item1 : solutions) {
            bool  isUnique = 1;
            for (auto item2 : all_solutions) {
                if (item1 == item2) {
                    isUnique = 0;
                    break;
                }
            }
            if (isUnique) {
                uniqueSols.push_back(item1);
            }
        }
        if (uniqueSols.size() == 0) {
            auto start = feasibilityCuts.cuts.rbegin();
            auto end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            auto start2 = optimalityCuts.cuts.rbegin();
            auto end2 = optimalityCuts.cuts.rend();
            for (; start2 != end2; ++start2) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
            }
            start = feasibilityCuts.cuts.rbegin();
            end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            start2 = optimalityCuts.cuts.rbegin();
            end2 = optimalityCuts.cuts.rend();
            for (; start2 != end2; ++start2) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
            }
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            CutCoefficients Y_bar_coef;
            double rhs = 0.0;
            CutType type = FEASIBILITY;
            Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
            auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
            if (feasibilityCuts.cuts.size()>1) {
                specialFeasCut  = feasibilityCuts.cuts[0];
                specialFeasCut.RHS = 0.0;
                for(auto item : specialFeasCut.cutCoeff) {
                    specialFeasCut.cutCoeff[item.first] = 0.0;
                }
                for (auto oneCut : feasibilityCuts.cuts) {
                    specialFeasCut.RHS += oneCut.RHS* 1.0 ;
                    for(auto item : oneCut.cutCoeff) {
                        specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
                    }
                }
                // cout << " app special feas: " ;
                if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
                    return {0,0,{}, false};
                }
            }

            /// opt Cuts
            CutCoefficients Y_bar_coef1;
            double rhs1 = 0.0;
            CutType type1 = OPTIMALITY ;
            Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
            auto TotaloptCutsNum = optimalityCuts.cuts.size();
            if (optimalityCuts.cuts.size()>1) {
                specialOptCut  = optimalityCuts.cuts[0];
                specialOptCut.RHS = 0.0;
                for(auto item : specialOptCut.cutCoeff) {
                    specialOptCut.cutCoeff[item.first] = 0.0;
                }
                for (auto oneCut : optimalityCuts.cuts) {
                    specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
                    for(auto item : oneCut.cutCoeff) {
                        specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
                    }
                }
                // cout << " app special opt: " ;
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
            }
            /// mix of special opt and special feas
            CutCoefficients Y_bar_coef2;
            double rhs2 = 0;
            CutType type2 = OPTIMALITY;
            Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
            if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
                verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
                for(auto item : specialFeasCut.cutCoeff) {
                    verySpecialCut.cutCoeff[item.first] = 0.0;
                }
                for(auto item : specialFeasCut.cutCoeff) {
                    verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
                }
                for(auto item : specialOptCut.cutCoeff) {
                    verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
                }
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
            }
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            /////////////////////////////////////////////////////
            start = feasibilityCuts.cuts.rbegin();
            end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            start2 = optimalityCuts.cuts.rbegin();
            end2 = optimalityCuts.cuts.rend();
            for (; start2 != end2; ++start2) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
            }
            start = feasibilityCuts.cuts.rbegin();
            end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            start2 = optimalityCuts.cuts.rbegin();
            end2 = optimalityCuts.cuts.rend();
            for (; start2 != end2; ++start2) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
            }

            auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
            for (auto& c_node: cutset) {
                c_node.ub = upperBound;
                c_node.lb = optimalLB;
            }
            return {lowerBound,upperBound,cutset, true};
        }

        for (auto solution : uniqueSols ) {
            all_solutions.push_back(solution);
        }
        // CutCoefficients Y_bar_coef;
        // double rhs = 0.0;
        // CutType type = FEASIBILITY;
        // Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
        // auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
        // if (feasibilityCuts.cuts.size()>1) {
        //     specialFeasCut  = feasibilityCuts.cuts[0];
        //     specialFeasCut.RHS = 0.0;
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         specialFeasCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for (auto oneCut : feasibilityCuts.cuts) {
        //         specialFeasCut.RHS += oneCut.RHS* 1.0 ;
        //         for(auto item : oneCut.cutCoeff) {
        //             specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
        //         }
        //     }
        //     // cout << " app special feas: " ;
        //     if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
        //         return {0,0,{}, false};
        //     }
        // }
        //
        // /// opt Cuts
        // CutCoefficients Y_bar_coef1;
        // double rhs1 = 0.0;
        // CutType type1 = OPTIMALITY ;
        // Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
        // auto TotaloptCutsNum = optimalityCuts.cuts.size();
        // if (optimalityCuts.cuts.size()>1) {
        //     specialOptCut  = optimalityCuts.cuts[0];
        //     specialOptCut.RHS = 0.0;
        //     for(auto item : specialOptCut.cutCoeff) {
        //         specialOptCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for (auto oneCut : optimalityCuts.cuts) {
        //         specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
        //         for(auto item : oneCut.cutCoeff) {
        //             specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
        //         }
        //     }
        //     // cout << " app special opt: " ;
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        // }
        // //////   mix of special opt and special feas
        // CutCoefficients Y_bar_coef2;
        // double rhs2 = 0;
        // CutType type2 = OPTIMALITY;
        // Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
        // if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
        //     verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] = 0.0;
        //     }
        //     for(auto item : specialFeasCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
        //     }
        //     for(auto item : specialOptCut.cutCoeff) {
        //         verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
        //     }
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        // }
        cout << " us:" << uniqueSols.size() << " ";
        double  newUpperBound;
        for (auto solution : uniqueSols ) {
            // cout << "got here! " ;
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
                if (!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
                // if (optimalityCuts.cuts.size()>0) {
                //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(lastEffectiveOptCut, optimalLB,upperBound),upperBound);
                //     if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
                // }
            }else {
                optimalityCuts.insertCut(cut);
                newUpperBound =  relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                newUpperBound >= upperBound;
                lastEffectiveOptCut = cut;
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}

            }
            // if (newUpperBound >= upperBound) {
            //     auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
            //     for (auto& c_node: cutset) {
            //         c_node.ub = upperBound;
            //         c_node.lb = optimalLB;
            //     }
            //     return {lowerBound,upperBound,cutset, true};
            // }else {
            //     upperBound = newUpperBound;
            // }

        }
    }
}



OutObject NodeExplorer::processX3(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";
    double lowerBound = std::numeric_limits<int64_t>::min();
    double upperBound = node.ub;
    DD relaxedDD{networkPtr, EXACT};
    // vector<vi> allSolutions;
    ///////////////////////////////////
    ///////////////////////////////////
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    relaxedDD.build(root2);
    if (relaxedDD.RelaxedisExact) {
        // cout << " exact " << endl;
        // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
        //     cout << i << " - " << relaxedDD.tree[i].size() << " global thing: " << relaxedDD.nodes[relaxedDD.tree[i][0]].globalLayer << endl;
        // }
        // cout << " 1 " ;
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            // cout << "zart" << endl;
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                cout << " pr.b.feas ";
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        // cout << " 2 " ;
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound = relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound);
            if (upperBound <= optimalLB) {
                cout << " pr.b.opt ";
                return {lowerBound,upperBound,{}, false};
            }
        }
        vector<vi> allSolutions={};
        while (true) {
            auto solution = relaxedDD.solution();
            if (std::find(allSolutions.begin(),allSolutions.end(),solution) != allSolutions.end()) {
                cout << "upperBound" << upperBound << endl;
                return {upperBound,upperBound,{}, true};
            }else {
                allSolutions.push_back(solution);
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            // cout << " 1 " << endl;
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
            }
            // cout << " 2 " << endl;
        }

    }else {
        // cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";
        // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
        //     cout << i << " - " << relaxedDD.tree[i].size() << endl;
        // }
        // cout << " d0 " << endl;
        // auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
        // cout << " before " ;
        // cout << cutset.size() << " ";
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        // cout << " d0 " << endl;
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        // cout << " d1 " << endl;
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }
        // start = feasibilityCuts.cuts.rbegin();
        // end = feasibilityCuts.cuts.rend();
        // for (; start != end; ++start) {
        //     if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
        //         return {std::numeric_limits<int64_t>::min(),0,{}, false};
        //     }
        // }
        // // cout << " d3 " << endl;
        // start2 = optimalityCuts.cuts.rbegin();
        // end2 = optimalityCuts.cuts.rend();
        // for (; start2 != end2; ++start2) {
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        // }
        // cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
        // cout << cutset.size() << " ";
        // cout << " trd " ;
        // start = feasibilityCuts.cuts.rbegin();
        // end = feasibilityCuts.cuts.rend();
        // for (; start != end; ++start) {
        //     if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
        //         return {std::numeric_limits<int64_t>::min(),0,{}, false};
        //     }
        // }
        // // cout << " d3 " << endl;
        // start2 = optimalityCuts.cuts.rbegin();
        // end2 = optimalityCuts.cuts.rend();
        // for (; start2 != end2; ++start2) {
        //     upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
        //     if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        // }
        // cout << " d4 " << endl;
        auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
        cout << cutset.size() << " ";
        for (auto& c_node: cutset) {
            c_node.ub = upperBound;
            c_node.lb = optimalLB;
        }
        return {lowerBound,upperBound,cutset, true};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if(optimalityCuts.cuts.size() == 0) {
            auto path = relaxedDD.solution();
            if (std::find(allSolutions.begin(), allSolutions.end(), path) == allSolutions.end()) {
                allSolutions.push_back(path);
                auto y_bar = w2y(path , networkPtr);
                GuroSolver solver{networkPtr, env};
                auto cut = solver.solveSubProblemInstance(y_bar);
                if (cut.cutType == FEASIBILITY) {
                    feasibilityCuts.insertCut(cut);
                    if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                        return {lowerBound,0,{}, false};
                    }
                }else {
                    optimalityCuts.insertCut(cut);
                    upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
                }

            }else {
                auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
                for (auto& c_node: cutset) {
                    c_node.ub = upperBound;
                    c_node.lb = optimalLB;
                }
                return {lowerBound,upperBound,cutset, true};
            }
        }else {
            auto start = feasibilityCuts.cuts.rbegin();
            auto end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            auto start2 = optimalityCuts.cuts.rbegin();
            auto end2 = optimalityCuts.cuts.rend();
            for (; start2 != end2; ++start2) {
                upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
            }
            start = feasibilityCuts.cuts.rbegin();
            end = feasibilityCuts.cuts.rend();
            for (; start != end; ++start) {
                if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                    return {std::numeric_limits<int64_t>::min(),0,{}, false};
                }
            }
            auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
            for (auto& c_node: cutset) {
                c_node.ub = upperBound;
                c_node.lb = optimalLB;
            }
            return {lowerBound,upperBound,cutset, true};
            while (true){
                set<vi> allPaths ={};
                vector<vi> uniquePaths ={};
                vi path;
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                ///////////////////////////////////////
                /// exponential number of solutions ///
                /// fix this part with something    ///
                ///////////////////////////////////////
                for (; start2 != end2; ++start2) {
                    auto newUpperBound =  relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound);
                    if (newUpperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                    if (newUpperBound == upperBound) {
                        path = relaxedDD.solution();
                        allPaths.insert(path);
                    }
                }
                if (allPaths.size() > 0) {
                    for (auto item: allPaths) {
                        if (std::find(allSolutions.begin(), allSolutions.end(), item) == allSolutions.end()) {
                            uniquePaths.push_back(item);
                        }
                    }
                }
                cout << "us:" << uniquePaths.size();
                if (uniquePaths.size() == 0) {
                    auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
                    for (auto& c_node: cutset) {
                        c_node.ub = upperBound;
                        c_node.lb = optimalLB;
                    }
                    return {lowerBound,upperBound,cutset, true};
                }
                for (auto item: uniquePaths) {
                    allSolutions.push_back(item);
                }
                /// we have the unique paths now
                for (auto item : uniquePaths) {
                    auto y_bar = w2y(item , networkPtr);
                    GuroSolver solver{networkPtr, env};
                    auto cut = solver.solveSubProblemInstance(y_bar);
                    if (cut.cutType == FEASIBILITY) {
                        feasibilityCuts.insertCut(cut);
                        if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                            return {lowerBound,0,{}, false};
                        }
                    }else {
                        optimalityCuts.insertCut(cut);
                        upperBound = min(relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound),upperBound);
                        if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
                    }
                }
                start = feasibilityCuts.cuts.rbegin();
                end = feasibilityCuts.cuts.rend();
                for (; start != end; ++start) {
                    if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                        return {std::numeric_limits<int64_t>::min(),0,{}, false};
                    }
                }
                start2 = optimalityCuts.cuts.rbegin();
                end2 = optimalityCuts.cuts.rend();
                for (; start2 != end2; ++start2) {
                    upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
                    if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
                }
                // auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
                // for (auto& c_node: cutset) {
                //     c_node.ub = upperBound;
                //     c_node.lb = optimalLB;
                // }
                // return {lowerBound,upperBound,cutset, true};
            }
        }
    }
}



OutObject NodeExplorer::processX4(const Node_t node, const double optimalLB) {
    cout << "layer: " << node.globalLayer << " ,ub: " << node.ub  << " ,cur opt: " << optimalLB << "    ";
    double lowerBound = std::numeric_limits<int64_t>::min();
    double upperBound = node.ub;
    DD relaxedDD{networkPtr, EXACT};
    DDNode root2{0, node.globalLayer, node.states, node.solutionVector};
    relaxedDD.build(root2);
    if (relaxedDD.RelaxedisExact) {
        // cout << " exact " << endl;
        // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
        //     cout << i << " - " << relaxedDD.tree[i].size() << " global thing: " << relaxedDD.nodes[relaxedDD.tree[i][0]].globalLayer << endl;
        // }
        // cout << " 1 " ;
        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            // cout << "zart" << endl;
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                cout << " pr.b.feas ";
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }
        // cout << " 2 " ;
        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound = relaxedDD.applyOptimalityCutHeuristic(*start2,optimalLB,upperBound);
            if (upperBound <= optimalLB) {
                cout << " pr.b.opt ";
                return {lowerBound,upperBound,{}, false};
            }
        }
        vector<vi> allSolutions={};
        while (true) {
            auto solution = relaxedDD.solution();
            if (std::find(allSolutions.begin(),allSolutions.end(),solution) != allSolutions.end()) {
                cout << "upperBound" << upperBound << endl;
                return {upperBound,upperBound,{}, true};
            }else {
                allSolutions.push_back(solution);
            }
            auto y_bar = w2y(solution , networkPtr);
            GuroSolver solver{networkPtr, env};
            auto cut = solver.solveSubProblemInstance(y_bar);
            // cout << " 1 " << endl;
            if (cut.cutType == FEASIBILITY) {
                feasibilityCuts.insertCut(cut);
                if(!relaxedDD.applyFeasibilityCutHeuristic(cut)) {
                    return {lowerBound,0,{}, false};
                }
            }else {
                optimalityCuts.insertCut(cut);
                upperBound = relaxedDD.applyOptimalityCutHeuristic(cut,optimalLB,upperBound);
                if (upperBound <= optimalLB) {return {lowerBound,0,{}, false};}
            }

        }
    }else {
        // cout << " exact " << endl;
        // for(int i = 0 ; i< relaxedDD.tree.size(); i++) {
        //     cout << i << " - " << relaxedDD.tree[i].size() << " global thing: " << relaxedDD.nodes[relaxedDD.tree[i][0]].globalLayer << endl;
        // }
        CutCoefficients Y_bar_coef;
        double rhs = 0.0;
        CutType type = FEASIBILITY;
        Cut specialFeasCut = Cut{type, rhs, Y_bar_coef};
        auto TotalFeasCutsNum = feasibilityCuts.cuts.size();
        if (feasibilityCuts.cuts.size()>1) {
            specialFeasCut  = feasibilityCuts.cuts[0];
            specialFeasCut.RHS = 0.0;
            for(auto item : specialFeasCut.cutCoeff) {
                specialFeasCut.cutCoeff[item.first] = 0.0;
            }
            for (auto oneCut : feasibilityCuts.cuts) {
                specialFeasCut.RHS += oneCut.RHS* 1.0 ;
                for(auto item : oneCut.cutCoeff) {
                    specialFeasCut.cutCoeff[item.first] += item.second * 1.0 ;
                }
            }
            // cout << " app special feas: " ;
            if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
                return {0,0,{}, false};
            }
        }

        /// opt Cuts
        CutCoefficients Y_bar_coef1;
        double rhs1 = 0.0;
        CutType type1 = OPTIMALITY ;
        Cut specialOptCut = Cut{type1, rhs1, Y_bar_coef1};
        auto TotaloptCutsNum = optimalityCuts.cuts.size();
        if (optimalityCuts.cuts.size()>1) {
            specialOptCut  = optimalityCuts.cuts[0];
            specialOptCut.RHS = 0.0;
            for(auto item : specialOptCut.cutCoeff) {
                specialOptCut.cutCoeff[item.first] = 0.0;
            }
            for (auto oneCut : optimalityCuts.cuts) {
                specialOptCut.RHS += oneCut.RHS * 1.0 /TotaloptCutsNum;
                for(auto item : oneCut.cutCoeff) {
                    specialOptCut.cutCoeff[item.first] += item.second * 1.0 / TotaloptCutsNum;
                }
            }
            // cout << " app special opt: " ;
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        }
        /// mix of special opt and special feas
        CutCoefficients Y_bar_coef2;
        double rhs2 = 0;
        CutType type2 = OPTIMALITY;
        Cut verySpecialCut = Cut{type1, rhs1, Y_bar_coef1};
        if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
            verySpecialCut.RHS  = specialOptCut.RHS + specialFeasCut.RHS;
            for(auto item : specialFeasCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] = 0.0;
            }
            for(auto item : specialFeasCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] += specialFeasCut.cutCoeff[item.first];
            }
            for(auto item : specialOptCut.cutCoeff) {
                verySpecialCut.cutCoeff[item.first] += specialOptCut.cutCoeff[item.first];
            }
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        }


        auto start = feasibilityCuts.cuts.rbegin();
        auto end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }

        auto start2 = optimalityCuts.cuts.rbegin();
        auto end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }

        if (optimalityCuts.cuts.size()>1 && feasibilityCuts.cuts.size()>1) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(specialFeasCut)) {
                return {0,0,{}, false};
            }

            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(specialOptCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}

            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(verySpecialCut, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {node.lb,upperBound,{}, false};}
        }

        cout << " sec tm " ;
        start = feasibilityCuts.cuts.rbegin();
        end = feasibilityCuts.cuts.rend();
        for (; start != end; ++start) {
            if(!relaxedDD.applyFeasibilityCutHeuristic(*start)) {
                return {std::numeric_limits<int64_t>::min(),0,{}, false};
            }
        }

        start2 = optimalityCuts.cuts.rbegin();
        end2 = optimalityCuts.cuts.rend();
        for (; start2 != end2; ++start2) {
            upperBound =  min(relaxedDD.applyOptimalityCutHeuristic(*start2, optimalLB,upperBound),upperBound);
            if (upperBound <= optimalLB) {return {lowerBound,upperBound,{}, false};}
        }

        auto cutset = relaxedDD.getExactCutsetRelaxed(upperBound);
        cout << cutset.size() << " ";
        for (auto& c_node: cutset) {
            c_node.ub = upperBound;
            c_node.lb = optimalLB;
        }
        return {lowerBound,upperBound,cutset, true};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}