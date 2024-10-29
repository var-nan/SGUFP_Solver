//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"

void DDSolver::NodeQueue::pushNodes(vector<Node_t> nodes) {
    // push bunch of nodes.
    q.insert(q.end(), nodes.begin(), nodes.end());
}

void DDSolver::NodeQueue::pushNode(Node_t&& node) {
    q.push_back(node);
}

Node_t DDSolver::NodeQueue::getNode() {
    auto node = q.back();
    q.pop_back();
    return node;
}

vector<Node_t> DDSolver::NodeQueue::getNodes(size_t n = 8) {
    vector<Node_t> nodes;

    for (size_t i = 0; i < n; i++) {
        nodes.push_back(q.back());
        q.pop_back();
    }
    return nodes;
}

Node_t DDSolver::getNode() {
    return nodeQueue.getNode();
}

void DDSolver::startSolve(const Network& network) {

    /*
     * 1. get node from node's queue
     * 2. If LB is better than node's UB, get rid of node and get another.
     * 4. Start NodeProcessor();
     * 5. get cutset and LB, and UB. update the respective global values.
     * 6. insert cutset nodes to the queue.
     */

    NodeExplorer explorer;

    while (!nodeQueue.empty()) { // conditional wait in parallel version

        Node_t node = nodeQueue.getNode();

        if (node.ub < getLB()) {
            #ifdef DEBUG
            cout << "selected node's upper bound < optimal lower bound." << endl;
            numNodesDiscarded++;
            #endif
            continue; // look for another
        }
        // start node processor
        auto result = explorer.process(network, node); // use co-routines to update globalLB in between.

        #ifdef DEBUG
        numNodesExplored++;
        numNodesUnnecessary += !result.success;
        #endif
        // either this node returns cutset or nothing.
        if (result.success) {
            if (result.lb > getLB()) {
                setLB(result.lb);
                nodeQueue.pushNodes(result.nodes);
            }
        }
    }

    #ifdef DEBUG
    assert(nodeQueue.empty());
    cout << "Work queue is empty." << endl;
    cout << "Optimal Solution: " << optimalLB << endl;
    displayStats();
    #endif
}

void DDSolver::start() {
    // startSolve();
}

double DDSolver::getLB() const{
    return optimalLB;
}

void DDSolver::setLB(double lb) {
    optimalLB = lb;
}