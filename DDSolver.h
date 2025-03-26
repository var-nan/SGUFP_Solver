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
#include "optimized.h"
#include <bit>
#include <gsl/gsl_statistics_double.h>

#ifndef POLL_FREQUENCY
#define POLL_FREQUENCY 512
#endif

#ifndef LOCAL_CUTS_LIMIT
#define LOCAL_CUTS_LIMIT 16
#endif

const unsigned int NUM_WORKERS = 3;
//const uint8_t SHIFT = 5;
// const uint16_t POLL_FREQUENC;

#define CUT_CONTAINER_CAPACITY 64
#define F_CUT_CACHE_SIZE 8
#define O_CUT_CACHE_SIZE 8
#define PROPORTION_OF_SHARE 0.4
#define PAYLOAD_CHECK_TICKS 10


//
// class Payload {
//
// public:
//     Payload() = default;
//     vector<Node_t> nodes_;
//     // std::atomic<uint8_t> status{}; // worker's current status.
//     volatile uint8_t status = 0;
//     std::mutex lock; // around nodes vector.
//     std::condition_variable cv; // to wake up the worker waiting for nodes.
//     std::atomic<uint> payloadStatus = 0;
//     CutContainer feasibilityCuts_;
//     CutContainer optimalityCuts_;
//
//     vector<Node_t> getNodes(bool &done); // called by worker.
//     void addNodesToWorker(vector<Node_t> nodes); // called by master.
//     bool masterRequireNodes() const noexcept; // called by master.
//     void askWorkerForNodes(); // called by master.
//     uint8_t getStatus() const noexcept;
//     vector<Node_t> getNodesFromWorker();
//     void addNodesToMaster(vector<Node_t> nodes);
//     void setStatus(uint8_t status_);
//
//      enum STATUS {
//         WORKER_NEEDS_NODES = 0x1,
//         MASTER_NEEDS_NODES = 0x2,
//         WORKER_WORKING = 0x4,
//         MASTER_ASSIGNED_NODES = 0x8,
//         WORKER_SHARED_NODES = 0x10,
//         NOT_ENOUGH_NODES_TO_SHARE = 0x20,
//         SOLVER_FINISHED = 0x40,
//         MASTER_RECEIVED_NODES = 0x80
//     };
//
// };
//
// class DDSolver {
//
//     class WorkerElement {
//         std::condition_variable cv;
//         std::mutex lock;
//         vector<Node_t> nodes_;
//         std::atomic<uint8_t> status{};
//         /*
//          *  0   : worker working in progress.
//          *  1   : worker needs nodes.
//          *  2   : master needs nodes.
//          *  4   : master assigned nodes to worker.
//          *  8   : worker assigned nodes to master.
//          *  16  : not enough nodes to return to master.
//          */
//
//     public:
//         WorkerElement() = default;
//
//         /**
//         * Worker calls this function to collect the nodes from master.
//         */
//         vector<Node_t> getWork() {
//
//             {
//                 std::scoped_lock l{lock};
//                 if (!nodes_.empty()) {
//                     /* either master placed some nodes initially, or worker placed some nodes previously
//                        for master on master's request. Status flag should be zero. */
//                     auto nodes = move(nodes_);
//                     // status.store(0, memory_order_release); // really necessary?
//                     return nodes;
//                 }
//             }
//             // indicate master and wait.
//             status.store(1,memory_order_release);
//             while (true) {
//                 std::unique_lock ul{lock};
//                 cv.wait(ul, [&]{return !nodes_.empty();});
//                 vector<Node_t> work = std::move(nodes_);
//                 status.store(0b0, memory_order_release);
//                 // nodes_.clear();
//                 return work;
//             }
//         }
//
//         /**
//         * Master adds the nodes to the element.
//         */
//         void addWork(vector<Node_t> work) {
//             // should be done by the master.
//             {
//                 std::unique_lock<std::mutex> ul{lock};
//                 nodes_ = move(work);
//             }
//             status.store(4, memory_order_release);
//             cv.notify_one();
//         }
//
//         bool submitWorkRequest() { // called by master.
//             // if worker itself require work?
//             auto val = status.load(memory_order_acquire);
//             if (val & 0b1000) return false;
//             status.store(2, memory_order_release);
//             return true;
//         }
//
//     };
//
//     class NodeQueue {
//         struct comparator {
//             bool operator() (const Node_t& node1, const Node_t& node2) const {
//                 // return (node1.ub + node2.lb) < (node2.ub + node2.lb); // need to optimize this with ints
//                 return node1.globalLayer > node2.globalLayer;
//             }
//         };
//         // use either queue or vector or priority queue.
//         // use mutex
//         priority_queue<Node_t, vector<Node_t>, comparator> q;
//         // stack<Node_t> q;
//     public:
//         NodeQueue() = default;
//         explicit NodeQueue(vector<Node_t> nodes): q{nodes.begin(), nodes.end()}{}
//
//         void pushNodes(vector<Node_t> nodes);
//         void pushNode(Node_t node);
//         Node_t getNode();
//         vector<Node_t> getNodes(size_t n);
//
//         [[nodiscard]] bool empty() const { return q.empty();}
//         [[nodiscard]] size_t size() const {return q.size();}
//
//     };
//
//     class Worker {
//         // some performance coutners
//         uint id;
//         vector<CutContainer *> oCutsGlobal;
//         vector<CutContainer *> fCutsGlobal;
//
//         constexpr static auto is_poll_time = [](const size_t processed) {
//             auto lsbs = processed >> std::bit_width(static_cast<uint8_t>(POLL_FREQUENCY));
//             auto res = lsbs ^ POLL_FREQUENCY;
//             return !res;
//         };
//
//         void shareCutsWithMaster(NodeExplorer& explorer, Payload& payload);
//
//     public:
//         explicit Worker(uint id_):id{id_}{};
//         void operator()(DDSolver& solver);
//     };
//
//     class Master {
//
//         NodeQueue nodeQueue;
//         vector<CutContainer *> oCutsGlobal;
//         vector<CutContainer *> fCutsGlobal;
//
//
//         void addCutsToGlobal(DDSolver &solver);
//
//     public:
//         Master() = default;
//         void operator()(DDSolver& solver);
//     };
//     NodeQueue nodeQueue; // global queue.
//     double optimalLB;
//
//     const shared_ptr<Network> networkPtr;
//
//     #ifdef SOLVER_STATS
//     size_t numNodesExplored = 0;
//     size_t numNodesFound = 0;
//     size_t numPrunedByBound = 0;
//     size_t numNodesUnnecessary = 0;
//     size_t numQueueEntered = 0;
//     void displayStats() const {
//         cout << "************************ Stats for nerds **************************" << endl;
//         cout << "Total number of nodes added to queue: " << numQueueEntered << endl;
//         cout << "Number of nodes processed: " << numNodesExplored << endl;
//         cout << "Number of nodes discarded by feasibility: " << numPrunedByBound << endl;
//         cout << "Number of unnecessarily processed nodes: " << numNodesUnnecessary << endl;
//         cout << "********************************************************************" << endl;
//     }
//     #endif
//
//     void process(NodeExplorer explorer);
// 	void processWork(unsigned int id, pair<CutContainer, CutContainer> cuts);
//     void processWork2(unsigned int id, pair<CutContainer, CutContainer> cuts);
//     void startMaster2();
//     void startMaster();
//     void startMaster3();
//     void processWork3(unsigned int id, pair<CutContainer, CutContainer> cuts);
//
// 	std::mutex queueLock;
// 	std::atomic<double> globalLB{numeric_limits<double>::lowest()};
//     std::atomic_bool isCompleted{false};
//     vector<Payload> workers{NUM_WORKERS};
//     Inavap::CutResource cutResources; //
//
//
//
// public:
//     explicit DDSolver(const shared_ptr<Network>& networkPtr_, const uint nWorkers):networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()}{}
//     explicit DDSolver(const shared_ptr<Network>& networkPtr_):networkPtr{networkPtr_}, optimalLB{std::numeric_limits<double>::lowest()} { }
//     // DDSolver() : optimalLB{std::numeric_limits<double>::lowest()}{}
//
//     [[nodiscard]] Node_t getNode();
//
//     [[nodiscard]] double getOptimalLB() const;
//     void setLB(double lb);
//
//     void initialize();
//     void start();
//     void startSolve(optional<pair<CutContainer, CutContainer>> initialCuts);
//     void startSolveParallel(optional<pair<CutContainer,CutContainer>> initialCuts);
//     pair<CutContainer, CutContainer> initializeCuts();
//     pair<CutContainer, CutContainer> initializeCuts2(size_t n = 50);
//
//     void startPThreadSolver();
// };


namespace Inavap {

    class DDSolver {

        friend class Worker;
        friend class Master;

        class Payload {
        public:
            enum STATUS {
                WORKER_WORKING = 0X1,
                WORKER_NEEDS_NODES = 0X2,
                WORKER_SHARED_NODES = 0X4,
                NOT_ENOUGH_NODES_TO_SHARE = 0X8,
                MASTER_NEEDS_NODES = 0X10,
                MASTER_ASSIGNED_NODES = 0X20,
                MASTER_RECEIVED_NODES = 0X40,
                SOLVER_FINISHED = 0X80
            };
            vector<Node> nodes;
            uint id = 0;
            std::mutex lock;
            std::condition_variable cv;
            std::atomic<uint> payloadStatus;
            CutContainer fCuts;
            CutContainer oCuts;

            Payload() = default;
            vector<Node> getNodes(uint8_t &done);
            void addCuts(optional<CutContainer> feasCuts, optional<CutContainer> optCuts);

        };

        class NodeQueue {

            /* To dynamically change the queue order at runtime, a new 8-bit field "priority" is added
             * to Node type. The priority value is calculated as:
             *
             * weight_DFS = 0.5, weight_BFS = 0.5
             *
             * priority = (weight_DFS) * globalLayer + (weight_BFS)*(1/globalLayer)
             *
             * The value is normalized to 2^8 if needed. The master thread can modify the priority by
             * changing the weights of DFS and BFS.
             */

            struct comparator {
                bool operator() (const Node &node1, const Node &node2) const {
                    return node1.globalLayer < node2.globalLayer; // DEPTH_FIRST_SEARCH
                }
            };

            priority_queue<Node, vector<Node>, comparator> pq;
        public:
            NodeQueue() = default;
            explicit NodeQueue(vector<Node> nodes) : pq{nodes.begin(), nodes.end()}{}

            void pushNode(Node node);
            void pushNodes(vector<Node> nodes);
            [[nodiscard]] Node popNode();
            [[nodiscard]] vector<Node> popNodes(size_t n);

            [[nodiscard]] bool empty() const {return pq.empty();}
            [[nodiscard]] size_t size() const {return pq.size();}
        };

        class Worker {
            const uint id;
        public:
            vector<CutContainer *> optCutsGlobal; // pointers to global optimality cut containers.
            vector<CutContainer *> feasCutsGlobal; // pointers to global feasibility cut containers.
            shared_ptr<Network> networkPtr;

            void shareCutsWithMaster(NodeExplorer& explorer, Inavap::DDSolver::Payload& payload);

        // public:
            explicit Worker(uint id_, const shared_ptr<Network>& networkPtr_): id{id_}, networkPtr{networkPtr_}{}
            void startWorker(DDSolver *solver);
            void printStats() const noexcept {
                auto str = "Worker: " + std::to_string(id) +" . FCC: " + std::to_string(feasCutsGlobal.size())
                    + ". OCC: " + std::to_string(optCutsGlobal.size()) +"\n";
                cout << str;
            }

            size_t nProcessed = 0; // number of nodes processed successfully.
            size_t nReceived = 0; // number of nodes received from the master.
            size_t nFalseProcessed = 0; // node is infeasible.
            size_t nFeasibilityPruned = 0;
            size_t nOptimalityPruned = 0;
        };

        class Master {
            NodeQueue nodeQueue;
            vector<Cut> tempOptCuts;
            vector<Cut> tempFeasCuts; // to be published global later.
            vector<CutContainer *> optCutsGlobal; // pointers to global optimality cut containers.
            vector<CutContainer *> feasCutsGlobal; // pointers to global feasibility cut containers.
            shared_ptr<Network> networkPtr;

            void addCutsToGlobal(DDSolver &solver, pair<vector<CutContainer>, vector<CutContainer>> &cuts);
            // bool isCutsShared(Payload &worker);
            size_t processNodes(DDSolver &solver, NodeExplorer &explorer, size_t n);

        public:
            explicit Master(const shared_ptr<Network>& networkPtr_) : networkPtr{networkPtr_}{};
            void startMaster (DDSolver &solver);
        };

        vector<Payload> payloads; // individual payloads for worker threads.
        CutResource CutResources;
        shared_ptr<Network> networkPtr;
        std::atomic_bool isCompleted{false};
        const uint16_t N_WORKERS;
        vector<Worker> workersGroup;
        vector<thread> workers; // worker threads.

        atomic<double> optimal = std::atomic<double>(std::numeric_limits<double>::lowest());


    public:
        explicit DDSolver(const shared_ptr<Network> &networkPtr_, uint16_t nWorkers):
            networkPtr{networkPtr_}, N_WORKERS{nWorkers}, payloads{nWorkers} {
            workers.reserve(N_WORKERS);
            workersGroup.reserve(N_WORKERS);
            cout << "Number of workers: " << N_WORKERS << endl;
        }


        void startSolver(double optimal);

        void printWorkerStats() const noexcept {
            vector<double> processed;
            processed.reserve(N_WORKERS);

            size_t n = processed.size() * 12 +32;
            n = std::max(n, static_cast<size_t>(72));
            std::string asterisk(n,'-');
            cout << asterisk <<endl;

            cout << "Processed: ";
            for (const auto& worker : workersGroup) {
                processed.push_back(static_cast<double>(worker.nProcessed));
                cout << worker.nProcessed << "  ";
            }
            cout << endl;
            double mean = gsl_stats_mean(processed.data(), 1, processed.size());
            double absdev = gsl_stats_absdev(processed.data(), 1, processed.size());
            double min_val= gsl_stats_min(processed.data(),1 , processed.size());
            double max_val = gsl_stats_max(processed.data(), 1, processed.size());

            double total = std::accumulate(processed.begin(), processed.end(), 0.0);

            cout << "Total: " << total << "\t Mean: " << mean << "\t Deviation: " << absdev
                        <<"\t Min: " << min_val <<"\t Max: " << max_val
                        << endl;
            cout << asterisk << endl;

            cout << endl;
            // print cuts information.
            for (const auto& worker : workersGroup) {
                cout << "(" << worker.feasCutsGlobal.size() << ", " << worker.optCutsGlobal.size() << ")" << "  ";
            }
            cout << endl;

            cout << asterisk << endl;

            for (const auto&  worker: workersGroup) {
                cout << "[" << worker.nFeasibilityPruned << " , " << worker.nOptimalityPruned << "]" << endl;
            }
            cout << asterisk << endl;
        }
    };
}


#endif //DDSOLVER_H
