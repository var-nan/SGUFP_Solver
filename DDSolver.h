//
// Created by nandgate on 10/24/2024.
//

#ifndef DDSOLVER_H
#define DDSOLVER_H

#include <thread>
#include <mutex>
#include <atomic>
#include "DD.h"
#include "grb.h"
#include "NodeExplorer.h"
#include <queue>
#include <stack>



class DDSolver {

    struct comparator {
        bool operator() (const Node_t& node1, const Node_t& node2) const {
            // return (node1.ub + node2.lb) < (node2.ub + node2.lb); // need to optimize this with ints
            return node1.globalLayer > node2.globalLayer;
        }
    };

    class NodeQueue {
        // use either queue or vector or priority queue.
        // use mutex
        priority_queue<Node_t, vector<Node_t>, comparator> q;
        // stack<Node_t> q;
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

    const shared_ptr<Network> networkPtr;

    #ifdef SOLVER_STATS
    size_t numNodesExplored = 0;
    size_t numNodesFound = 0;
    size_t numPrunedByBound = 0;
    size_t numNodesUnnecessary = 0;
    size_t numQueueEntered = 0;
    void displayStats() const {
        cout << "************************ Stats for nerds **************************" << endl;
        cout << "Total number of nodes added to queue: " << numQueueEntered << endl;
        cout << "Number of nodes processed: " << numNodesExplored << endl;
        cout << "Number of nodes discarded by feasibility: " << numPrunedByBound << endl;
        cout << "Number of unnecessarily processed nodes: " << numNodesUnnecessary << endl;
        cout << "********************************************************************" << endl;
    }
    #endif

    void process(NodeExplorer explorer);
	void processWork();

	std::mutex queueLock;
	std::atomic<double> globalLB{numeric_limits<double>::lowest()};

public:


    explicit DDSolver(const shared_ptr<Network>& networkPtr_):networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()} { }
    // DDSolver() : optimalLB{std::numeric_limits<double>::lowest()}{}

    [[nodiscard]] Node_t getNode();

    [[nodiscard]] double getOptimalLB() const;
    void setLB(double lb);

    void initialize();
    void start();
    void startSolve(optional<pair<CutContainer, CutContainer>> initialCuts);
    pair<CutContainer, CutContainer> initializeCuts();
    pair<CutContainer, CutContainer> initializeCuts2(size_t n = 50);

	void startPThreadSolver();
};



#endif //DDSOLVER_H
