//
// Created by nandgate on 10/24/2024.
//

#ifndef DDSOLVER_H
#define DDSOLVER_H

#include "DD.h"
#include "grb.h"
#include "NodeExplorer.h"
#include <queue>




class DDSolver {

    struct comparator {
        bool operator() (const Node_t& node1, const Node_t& node2) const {
            return node1.ub > node2.ub;
        }
    };

    class NodeQueue {
        // use either queue or vector or priority queue.
        // use mutex
        priority_queue<Node_t, vector<Node_t>, comparator> q;
    public:
        NodeQueue() = default;

        void pushNodes(vector<Node_t> nodes);
        void pushNode(Node_t node);
        Node_t getNode();
        vector<Node_t> getNodes(size_t n);

        [[nodiscard]] bool empty() const { return q.empty();}
        [[nodiscard]] size_t size() const {return q.size();}

    };

    NodeQueue nodeQueue;
    double optimalLB;

    #ifdef DEBUG
    size_t numNodesExplored = 0;
    size_t numNodesFound = 0;
    size_t numPrunedByBound = 0;
    size_t numNodesUnnecessary = 0;
    void displayStats() const {
        cout << "************************ Stats for nerds **************************" << endl;
        cout << "Total number of nodes found: " << numNodesFound << endl;
        cout << "Number of nodes processed: " << numNodesExplored << endl;
        cout << "Number of nodes discarded by feasibility: " << numPrunedByBound << endl;
        cout << "Number of unnecessarily processed nodes: " << numNodesUnnecessary << endl;
        cout << "********************************************************************" << endl;
    }
    #endif


public:

    DDSolver():optimalLB{std::numeric_limits<double>::lowest()} { }
    [[nodiscard]] Node_t getNode();

    [[nodiscard]] double getOptimalLB() const;
    void setLB(double lb);

    void initialize();
    void start();
    void startSolve(const Network& network);
};



#endif //DDSOLVER_H
