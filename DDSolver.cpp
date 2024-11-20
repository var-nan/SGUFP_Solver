//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"
#include <random>
#include <chrono>
#include <omp.h>

void DDSolver::NodeQueue::pushNodes(vector<Node_t> nodes) {
    // push bunch of nodes.
    // q.insert(q.end(), nodes.begin(), nodes.end());
    for (auto& node : nodes) {
        q.push(node);
    }
}

void DDSolver::NodeQueue::pushNode(Node_t node) {
    q.push(node);
}

Node_t DDSolver::NodeQueue::getNode() {
    // auto node = q.back(); q.pop_back();
    // auto node = q.front(); q.pop();
    auto node = q.top(); q.pop();
    return node;
}

vector<Node_t> DDSolver::NodeQueue::getNodes(size_t n = 8) {
    vector<Node_t> nodes;
    //
    // for (size_t i = 0; i < n; i++) {
    //     nodes.push_back(q.back());
    //     q.pop_back();
    // }
    return nodes;
}

Node_t DDSolver::getNode() {
    return nodeQueue.getNode();
}

void DDSolver::initialize() {
    // place root node to queue
    Node_t node;

    node.lb = std::numeric_limits<int64_t>::min();
    node.ub = std::numeric_limits<int64_t>::max();
    node.globalLayer = 0;

    nodeQueue.pushNode(node);
}

void DDSolver::startSolve(optional<pair<CutContainer, CutContainer>> initialCuts = nullopt) {
    if (initialCuts) {
        NodeExplorer explorer{networkPtr, initialCuts.value()};
        process(explorer);
    }
    else {
        NodeExplorer explorer{networkPtr};
        process(explorer);
    }
}

void DDSolver::process(NodeExplorer explorer) {

    /*
     * 1. get node from node's queue
     * 2. If LB is better than node's UB, get rid of node and get another.
     * 4. Start NodeProcessor();
     * 5. get cutset and LB, and UB. update the respective global values.
     * 6. insert cutset nodes to the queue.
     */

    // NodeExplorer explorer{networkPtr, initialCuts.value()};

    while (!nodeQueue.empty()) { // conditional wait in parallel version

        Node_t node = nodeQueue.getNode();
        // cout << "Procesisng Node from layer: "<< node.globalLayer << " LB: " << node.lb << " , UB: " << node.ub << " global: " << getOptimalLB()<< endl;
        #ifdef DEBUG
        // cout << "Processing node from layer: " << node.globalLayer << " lb: " << node.lb << " , ub: " << node.ub;
        cout << " . global lower bound: " << getOptimalLB() << endl;
        #endif
        if (node.ub < getOptimalLB()) {
            #ifdef SOLVER_STATS
            numPrunedByBound++;
            // cout << "Pruned by bound." << endl;
            #endif
            continue; // look for another
        }

        // start node processor
        auto result = explorer.process3(node, getOptimalLB()); // use co-routines to update globalLB in between.

        #ifdef SOLVER_STATS
        numNodesExplored++;
        numNodesUnnecessary += !result.success;
        #endif
        // either this node returns cutset or nothing.
        if (result.success) {
            if (result.lb > getOptimalLB()) {
                setLB(result.lb);
                const auto now = std::chrono::system_clock::now();
                const auto t_c = std::chrono::system_clock::to_time_t(now);
                cout << "Optimal LB: " << getOptimalLB() <<" set at "<< std::ctime(&t_c) << endl;
            }
            if (!result.nodes.empty()) {
                if (result.ub > getOptimalLB()) {
                    nodeQueue.pushNodes(result.nodes);
                    #ifdef SOLVER_STATS
                    // cout << result.nodes.size() << " nodes entered queue." << endl;
                    numQueueEntered += result.nodes.size();
                    #endif
                }
                else {
#ifdef DEBUG
                    cout << "Upper bound is worse than the global lb: " << getOptimalLB() << endl;
#endif
                }
            }
        }
        // explorer.clearCuts();
    }

    #ifdef SOLVER_STATS
    displayStats();
    explorer.displayCutStats();
    #endif
    cout << "Optimal Solution: " << getOptimalLB() << endl;
}

void DDSolver::start() {
    // startSolve();
}

double DDSolver::getOptimalLB() const{
    double lb = 0;
    #pragma omp atomic read
	lb = optimalLB;

    return lb;
}

void DDSolver::setLB(double lb) {
    #pragma omp critical
	{
		optimalLB = (lb > optimalLB)? lb : optimalLB;
		#pragma omp flush(optimalLB)
    }
}

DDNode node2DDdfsNode(Node_t node) {
    DDNode newNode;
    newNode.states = unordered_set<int>(node.states.begin(), node.states.end());
    newNode.solutionVector = node.solutionVector;
    newNode.globalLayer = node.globalLayer;
    newNode.nodeLayer = 0;
    return newNode;
}

pair<CutContainer, CutContainer> DDSolver::initializeCuts2(size_t n) {

    DD dd{networkPtr, EXACT};
    DDNode root{0};
    dd.build(root);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> dist(0, std::numeric_limits<int>::max());

    // srand(time(nullptr));

    // start with root and find random paths to terminal.
    vector<vi> solutions;
    CutContainer fCuts{FEASIBILITY};
    CutContainer oCuts{OPTIMALITY};
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    solutions.reserve(n);

    while (n) {
        // get path from root to terminal.
        vi solution;
        ulint currentId = 0;

        for (size_t i = 0; i < dd.tree.size() -2; i++) {
            const auto& node = dd.nodes.at(currentId);
            const auto& nArcs = node.outgoingArcs.size();
            auto selection = node.outgoingArcs[ dist(gen)%nArcs];
            const auto& arc = dd.arcs.at(selection);
            solution.push_back(arc.decision);
            currentId = arc.head;
        }

        // check if solution exists.
        bool isExist = false;
        for (const auto& sol : solutions) {
            if (sol == solution) isExist = true; break;
        }
        if (isExist) continue;
        cout << "Solution selected: "; for (auto s: solution) cout << s <<" "; cout <<endl;
        solutions.push_back(solution);
        // cout << "Solution: "; for (auto sol: solution) cout << sol << " "; cout << endl;
        // get cut
        GuroSolver solver{networkPtr, env};
        const auto& y_bar = w2y(solution, networkPtr);
        auto cut = solver.solveSubProblemInstance(y_bar, 0); // LATER set random scenario.
        if (cut.cutType == FEASIBILITY) {
            fCuts.insertCut(cut);
        }
        else oCuts.insertCut(cut);
        n--;
    }
    return make_pair(fCuts, oCuts);
}

pair<CutContainer, CutContainer> DDSolver::initializeCuts() {
    // build a restricted tree, get exact cutsets, build subtrees and get solutions for each of subtree.

    // change the max_width
    constexpr size_t max_width = MAX_WIDTH;
    // #undef MAX_WIDTH
    // #define MAX_WIDTH 20
    // #undef DEBUG

    // cout << "Current max width is : " << MAX_WIDTH << endl;

    /// build tree
    DDNode root{0};
    root.nodeLayer = 0;
    root.globalLayer = 0;
    DD relaxed{networkPtr, EXACT};
    DDNode newRoot{0};
    relaxed.build(newRoot);
    DD restricted {networkPtr,RESTRICTED};
    auto cutset = restricted.build(root);
    cout << "Cutset size " << cutset.value().size() << endl;

    CutContainer fCuts{FEASIBILITY};
    CutContainer oCuts{OPTIMALITY};

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag,0);
    if (cutset) {

        for (const auto& node: cutset.value()) {
            DDNode ddNode = node2DDNode(node);
            DD subTree{networkPtr,RESTRICTED};
            auto _ = subTree.build(ddNode);

            // get solution
            auto sol = subTree.solution();
            // solve sub problem
            GuroSolver solver{networkPtr,env};
            auto y_bar = w2y(sol, networkPtr);
            auto cut = solver.solveSubProblemInstance(y_bar, 0);
            // add cut to containers
            if (cut.cutType == FEASIBILITY) {
                if (!fCuts.isCutExists(cut))
                    fCuts.insertCut(cut);
            }
            else {
                if(!oCuts.isCutExists(cut)) oCuts.insertCut(cut);
            }
        }

    }

    // #undef MAX_WIDTH
    // #define DEBUG
    // #define MAX_WIDTH max_width

    // for (auto cut : fCuts.cuts) {
    //     cout << endl << "RHS: " << cut.RHS << endl;
    //     if (cut.cutType == FEASIBILITY) cout << "Type: FEASIBILITY" << endl;
    //     else cout << "Type: OPTIMALITY" << endl;
    //     for (auto [k,v] : cut.cutCoeff) {
    //         auto [i,q,j] = k;
    //         cout << "i: " << i << ", q: " << q << " j: " << j << " val: " << v << endl;
    //     }
    // }

    // cout << "Current max width is : " << MAX_WIDTH  << endl;
    return make_pair(fCuts, oCuts);

}

pair<CutContainer, CutContainer> DDSolver::initializeCutsParallel(size_t n) {

    vector<vi> solutions;
    int temp_n = n;

    // vector<GRBEnv> environments; for (size_t i = 0; i < 4; i++) environments.emplace_back(GRBEnv());

    CutContainer fCuts{FEASIBILITY};
    CutContainer oCuts{OPTIMALITY};

    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dist(0, std::numeric_limits<int>::max());

        DDNode root{0};
        DD dd{networkPtr, EXACT};
        dd.build(root);

        while (n) {
            vi solution;
            ulint currentId = 0;

            for (size_t i = 0; i < dd.tree.size() -2; i++) {
                const auto& node = dd.nodes.at(currentId);
                const auto& nArcs = node.outgoingArcs.size();
                auto selection = node.outgoingArcs[ dist(gen)%nArcs];
                const auto& arc = dd.arcs.at(selection);
                solution.push_back(arc.decision);
                currentId = arc.head;
            }
            // check if solution exists.
            bool isExist = false;
            for (const auto& sol : solutions) {
                if (sol == solution) {isExist = true; break;}
            }
            if (isExist) continue;
            solutions.push_back(solution);
            n--;
        }
    }
    n = temp_n;
    // generate cuts in parallel
    #pragma omp parallel num_threads(4) shared(fCuts, oCuts, solutions, n)
    {
        // cout << "Number of threads in parallel: " << omp_get_num_threads() << endl;
        // divide work between threads.
        int n_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();
        int chunk = n/n_threads;
        int start = chunk*thread_id;
        int end = min(start + chunk, static_cast<int>(n));

        CutContainer localFCuts{FEASIBILITY};
        CutContainer localOCuts{OPTIMALITY};

        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 0);

        for (int i = start; i < end; i++) {
            const auto& solution = solutions[i];
            auto y_bar = w2y(solution, networkPtr);
            GuroSolver solver{networkPtr, env};

            auto cut = solver.solveSubProblemInstance(y_bar, 0);

            if (cut.cutType == FEASIBILITY) {
                localFCuts.insertCut(cut);
            }
            else localOCuts.insertCut(cut);
        }

        #pragma omp critical
        {
            for (auto& cut : localFCuts.cuts) fCuts.insertCut(cut);
            for (auto& cut : localOCuts.cuts) oCuts.insertCut(cut);
        }
    }
    return {fCuts, oCuts};
}


void DDSolver::startSolveParallel(optional<pair<CutContainer, CutContainer>> initialCuts) {


    // create restricted tree and get cutset
    {
        DDNode root{0};
        DD dd{networkPtr, RESTRICTED};
        auto cutset = dd.build(root);
        // cutset should be valid.
        assert(cutset);
        nodeQueue.pushNodes(cutset.value());
    }

    vector<tuple<int,int,int>> thread_stats;

    int NUM_THREADS = 4;

    // distribute work by master.
    vector<NodeQueue> queues(NUM_THREADS);
    for (int i = 0; i < nodeQueue.size() ;i++) {
        queues[i%4].pushNode(nodeQueue.getNode());
    }

    #pragma omp parallel num_threads(NUM_THREADS) default(none) shared(queues, cout) firstprivate(initialCuts)
    {

        #pragma omp single
        {
            cout << "# threads: " << omp_get_num_threads() << endl;
        }

        int thread_id = omp_get_thread_num();
       // create local queue
        NodeQueue localQ; // assign comparing strategy later.
        while (!queues[thread_id].empty()) localQ.pushNode(queues[thread_id].getNode());
        int n_threads = omp_get_num_threads();

        #pragma omp barrier

        uint n_processed =0;
        uint n_pruned_by_bound = 0;

        NodeExplorer explorer{networkPtr, initialCuts.value()};
        // until work queue is empty.

        while (true) {
//			double optLB = 0;
            {
                Node_t node;
                bool isEmpty = false;
                #pragma omp critical
                {
                    if (nodeQueue.empty()) isEmpty = true;
                    else node = nodeQueue.getNode();
					#pragma omp flush (nodeQueue)
					// call getOptimalLB here.
                }
                if (localQ.empty() && isEmpty) break;
                localQ.pushNode(node);
            }

			double optLB = getOptimalLB(); // get latest value if processed a node.

            while (!localQ.empty()) {
                Node_t node = localQ.getNode();
//                double optLB = getOptimalLB();
                if (node.ub < optLB) {n_pruned_by_bound++; continue;}

                auto result = explorer.process3(node, optLB);
                n_processed++;

                if (result.success) {
                    optLB = getOptimalLB();
                    if (result.lb > optLB) {
                        // update opt
                        setLB(result.lb);
                        const auto now = std::chrono::system_clock::now();
                        const auto t_c = std::chrono::system_clock::to_time_t(now);
                        cout << "Thread : " << thread_id << ", Optimal LB: " << getOptimalLB() << " set at "
                            << std::ctime(&t_c)<< endl;
                    }
                    if(!result.nodes.empty()) {
//                        optLB = getOptimalLB();
                        if (result.ub > optLB) {
                            const int LOCAL_QUEUE_LIMIT = 4000;
                            auto result_nodes = result.nodes;
                            while ((localQ.size() < LOCAL_QUEUE_LIMIT) && !result_nodes.empty()) {
                                localQ.pushNode(result_nodes.back());
                                result_nodes.pop_back();
                            }
                            if (!result_nodes.empty()) {
								#pragma omp critical
	                            {
									nodeQueue.pushNodes(result_nodes);
									#pragma omp flush(nodeQueue)
								}
							}
                        }
                    }
                }
            }
        }

        // thread stats
        cout << "Thread " << thread_id << " , processed " << n_processed << " nodes." << endl;
        explorer.displayCutStats();
    }

    cout << "Optimal solution: " << getOptimalLB() << endl;

}