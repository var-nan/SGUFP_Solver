//
// Created by nandgate on 10/24/2024.
//

#ifndef DDSOLVER_H
#define DDSOLVER_H

#include <thread>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include "DD.h"
#include "grb.h"
#include "NodeExplorer.h"
#include <queue>
#include <stack>

class Payload {
    enum STATUS {
        WORKER_NEEDS_NODES = 0x1,
        MASTER_NEEDS_NODES = 0x2,
        WORKER_WORKING = 0x4,
        MASTER_ASSIGNED_NODES = 0x8,
        WORKER_SHARED_NODES = 0x100,
        NOT_ENOUGH_NODES_TO_SHARE = 0x10,
    };

    std::condition_variable cv;
    std::mutex lock;
    vector<Node_t> nodes_;
    std::atomic<uint8_t> status{};
    /*
     *  0   : worker working in progress.
     *  1   : worker needs nodes.
     *  2   : master needs nodes.
     *  4   : worker working in progress.
     *  8   : master assigned nodes to worker.
     *  16  : worker shared nodes to master (load balance).
     *  32  : not enough nodes to share to master.
     */

public:
    Payload() = default;

    vector<Node_t> getNodes();
    void addNodes(vector<Node_t> nodes);
    bool masterRequireNodes() const noexcept;


};

class DDSolver {

    friend class Worker;
    friend class Master;

    enum STATUS_FLAGS {
        ASSIGNED,
        REQUESTED = 8,
        REQUEST_MASTER = 1,
        REQUEST_WORKER = 2,
        PROCESSING = 0
    };

    class WorkerElement {
        std::condition_variable cv;
        std::mutex lock;
        vector<Node_t> nodes_;
        std::atomic<uint8_t> status{};
        /*
         *  0   : worker working in progress.
         *  1   : worker needs nodes.
         *  2   : master needs nodes.
         *  4   : master assigned nodes to worker.
         *  8   : worker assigned nodes to master.
         *  16  : not enough nodes to return to master.
         */

    public:
        WorkerElement() = default;

        /**
        * Worker calls this function to collect the nodes from master.
        */
        vector<Node_t> getWork() {

            {
                std::scoped_lock l{lock};
                if (!nodes_.empty()) {
                    /* either master placed some nodes initially, or worker placed some nodes previously
                       for master on master's request. Status flag should be zero. */
                    auto nodes = move(nodes_);
                    // status.store(0, memory_order_release); // really necessary?
                    return nodes;
                }
            }
            // indicate master and wait.
            status.store(1,memory_order_release);
            while (true) {
                std::unique_lock ul{lock};
                cv.wait(ul, [&]{return !nodes_.empty();});
                vector<Node_t> work = std::move(nodes_);
                status.store(0b0, memory_order_release);
                // nodes_.clear();
                return work;
            }
        }

        /**
        * Master adds the nodes to the element.
        */
        void addWork(vector<Node_t> work) {
            // should be done by the master.
            {
                std::unique_lock<std::mutex> ul{lock};
                nodes_ = move(work);
            }
            status.store(4, memory_order_release);
            cv.notify_one();
        }

        bool submitWorkRequest() { // called by master.
            // if worker itself require work?
            auto val = status.load(memory_order_acquire);
            if (val & 0b1000) return false;
            status.store(2, memory_order_release);
            return true;
        }

    };

    class NodeQueue {
        struct comparator {
            bool operator() (const Node_t& node1, const Node_t& node2) const {
                // return (node1.ub + node2.lb) < (node2.ub + node2.lb); // need to optimize this with ints
                // return 10* (50 - node1.globalLayer) + node1.ub > 10* (50-node2.globalLayer) + node2.ub; // birth first
                // return node1.ub -  node1.globalLayer  < node2.ub - node2.globalLayer ;
                // return node1.ub  < node2.ub  ;
                // return node1.globalLayer > node2.globalLayer; // birth first
                return node1.globalLayer < node2.globalLayer; // depth first
                // return node1.ub + 200 * node1.globalLayer < node2.ub + 200 * node2.globalLayer;
            }
        };
        // use either queue or vector or priority queue.
        // use mutex
        priority_queue<Node_t, vector<Node_t>, comparator> q;
        // stack<Node_t> q;
    public:
        NodeQueue() = default;
        explicit NodeQueue(vector<Node_t> nodes): q{nodes.begin(), nodes.end()}{}

        void pushNodes(vector<Node_t> nodes);
        void pushNode(Node_t node);
        Node_t getNode();
        vector<Node_t> getNodes(size_t n);

        [[nodiscard]] bool empty() const { return q.empty();}
        [[nodiscard]] size_t size() const {return q.size();}

    };

    NodeQueue nodeQueue; // global queue.
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
	void processWork(NodeQueue q, pair<CutContainer, CutContainer> cuts);

	std::mutex queueLock;
	std::atomic<double> globalLB{numeric_limits<double>::lowest()};
    std::atomic_bool isCompleted{false};
    vector<Payload> workers{2};


public:


    explicit DDSolver(const shared_ptr<Network>& networkPtr_):networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()} { }
    // DDSolver() : optimalLB{std::numeric_limits<double>::lowest()}{}

    [[nodiscard]] Node_t getNode();

    [[nodiscard]] double getOptimalLB() const;
    void setLB(double lb);

    void initialize(double opt);
    void start();
    void startSolve(optional<pair<CutContainer, CutContainer>> initialCuts);
    pair<CutContainer, CutContainer> initializeCuts();
    pair<CutContainer, CutContainer> initializeCuts2(size_t n = 50);
    pair<CutContainer, CutContainer> initializeCuts3();
    pair<CutContainer, CutContainer> initializeCuts4();
	void startPThreadSolver();
};



#endif //DDSOLVER_H
