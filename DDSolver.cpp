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

void DDSolver::startSolve() {

    /*
     * 1. get node from node's queue
     * 2. communicate LB with solver's current LB
     * 3. If LB is better than node's UB, get rid of node and get another.
     * 4. Start NodeProcessor();
     * 5. Update UB
     *
     */

    while (!nodeQueue.empty()) { // conditional wait in parallel version

        Node_t node = nodeQueue.getNode();

        if (node.ub < lowerBoundOptimal) {
            continue; // look for another
        }

        // start node processor
        NodeExplorer explorer;
        explorer.process(node);
    }


}

void DDSolver::start() {
    startSolve();
}

