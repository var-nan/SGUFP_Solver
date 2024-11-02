//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"

#include <chrono>

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
    auto node = q.front(); q.pop();
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

    node.lb = std::numeric_limits<double>::lowest();
    node.ub = std::numeric_limits<double>::max();
    node.globalLayer = 0;

    nodeQueue.pushNode(node);
}

void DDSolver::startSolve(optional<pair<CutContainer, CutContainer>> initialCuts = nullopt) {

    /*
     * 1. get node from node's queue
     * 2. If LB is better than node's UB, get rid of node and get another.
     * 4. Start NodeProcessor();
     * 5. get cutset and LB, and UB. update the respective global values.
     * 6. insert cutset nodes to the queue.
     */

    NodeExplorer explorer{networkPtr, initialCuts.value()};

    while (!nodeQueue.empty()) { // conditional wait in parallel version

        Node_t node = nodeQueue.getNode();
        // if (getOptimalLB() > 5700) break;

        if (node.ub < getOptimalLB()) {
            continue; // look for another
        }

        // start node processor
        auto result = explorer.process(node); // use co-routines to update globalLB in between.

        #ifdef DEBUG
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
                cout << "Global Lower bound so far: " << getOptimalLB() << endl;
            }
            if (!result.nodes.empty()) {
                if (result.ub > getOptimalLB()) {
                    nodeQueue.pushNodes(result.nodes);
                    #ifdef DEBUG
                    numQueueEntered += result.nodes.size();
                    #endif
                }
            }
        }
        // explorer.clearCuts();
    }

    #ifdef DEBUG
    //assert(nodeQueue.empty());
    cout << "Work queue is empty." << endl;
    cout << "Optimal Solution: " << getOptimalLB() << endl;
    displayStats();
    #endif
    cout << "Optimal Solution: " << getOptimalLB() << endl;
}

void DDSolver::start() {
    // startSolve();
}

double DDSolver::getOptimalLB() const{
    return optimalLB;
}

void DDSolver::setLB(double lb) {
    optimalLB = lb;
}

DDNode node2DDdfsNode(Node_t node) {
    DDNode newNode;
    newNode.states = unordered_set<int>(node.states.begin(), node.states.end());
    newNode.solutionVector = node.solutionVector;
    newNode.globalLayer = node.globalLayer;
    newNode.nodeLayer = 0;
    return newNode;
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
