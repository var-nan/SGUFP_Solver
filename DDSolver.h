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

static unsigned int NUM_WORKERS = 3;

enum STATUS {
        WORKER_NEEDS_NODES = 0x1,
        MASTER_NEEDS_NODES = 0x2,
        WORKER_WORKING = 0x4,
        MASTER_ASSIGNED_NODES = 0x8,
        WORKER_SHARED_NODES = 0x10,
        NOT_ENOUGH_NODES_TO_SHARE = 0x20,
        SOLVER_FINISHED = 0x40,
        MASTER_RECEIVED_NODES = 0x80
    };

class Payload {

    // std::condition_variable cv;
    // std::mutex lock;

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
    vector<Node_t> nodes_;
    // std::atomic<uint8_t> status{}; // worker's current status.
    volatile uint8_t status = 0;
    std::mutex lock; // around nodes vector.
    std::condition_variable cv; // to wake up the worker waiting for nodes.
    // uint8_t payloadStatus = 0;

    vector<Node_t> getNodes(bool &done); // called by worker.
    void addNodesToWorker(vector<Node_t> nodes); // called by master.
    bool masterRequireNodes() const noexcept; // called by master.
    void askWorkerForNodes(); // called by master.
    uint8_t getStatus() const noexcept;
    vector<Node_t> getNodesFromWorker();
    void addNodesToMaster(vector<Node_t> nodes);
    void setStatus(uint8_t status_);


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
                return node1.globalLayer > node2.globalLayer;
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
	void processWork(unsigned int id, pair<CutContainer, CutContainer> cuts);
    void processWork2(unsigned int id, pair<CutContainer, CutContainer> cuts);
    void startMaster2();
    void startMaster();
    void startMaster3();
    void processWork3(unsigned int id, pair<CutContainer, CutContainer> cuts);

    OutObject callProcess(NodeExplorer& explorer, Node_t& node, double opt, int process = 4);

	std::mutex queueLock;
	std::atomic<double> globalLB{numeric_limits<double>::lowest()};
    std::atomic_bool isCompleted{false};
    vector<Payload> workers;
    int process_number = 4;
    bool saveCuts = true;


public:


    explicit DDSolver(const shared_ptr<Network>& networkPtr_, int process_ = 4)
            :networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()}, process_number(process_) { }
    explicit DDSolver(const shared_ptr<Network>& networkPtr_, int process_ = 4, bool saveCuts_ = true)
            :networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()}, process_number(process_) , saveCuts{saveCuts_} {}
    // DDSolver() : optimalLB{std::numeric_limits<double>::lowest()}{}

    [[nodiscard]] Node_t getNode();

    [[nodiscard]] double getOptimalLB() const;
    void setLB(double lb);

    void initialize();
    void start();
    void startSolve(optional<pair<CutContainer, CutContainer>> initialCuts);
    pair<CutContainer, CutContainer> initializeCuts();
    pair<CutContainer, CutContainer> initializeCuts2(size_t n = 50);

	void startPThreadSolver(size_t nthreads = 4);
};



#endif //DDSOLVER_H
