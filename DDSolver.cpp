//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"
#include <random>
#include <chrono>
#include <omp.h>
#define OMP_NUM_THREADS 4
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

	while (!q.empty() && n) {
		//
		if (n--) {
			nodes.push_back(q.top()); q.pop();
		}
		// else break;
	}
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
    newNode.states = set<int>(node.states.begin(), node.states.end());
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


//void DDSolver::processWork() {
//	NodeExplorer explorer{networkPtr};
//
//	NodeQueue localQueue;
//
//	double zOpt = globalLB.load();
//	size_t nProcessed = 0;
//	Node_t rootNode;
//	rootNode.lb = numeric_limits<double>::lowest();
//	rootNode.ub = numeric_limits<double>::max();
//	rootNode.globalLayer = 0;
//	localQueue.pushNode(rootNode);
//	// insert all but one node to the local queue.
////	while (!nodeQueue.empty()) {localQueue.pushNode(nodeQueue.getNode());}
////	localQueue = nodeQueue;
//	while (true) {
//		// acquire lock get a node from global queue
////		Node_t node;
////		{
////			lock_guard<mutex> guard(queueLock);
////			if (nodeQueue.empty()) break;
////			node = nodeQueue.getNode();
////		}
////
////		localQueue.pushNode(node);
//		if (localQueue.empty()) break;
//		zOpt = globalLB.load(memory_order_acquire);
//
//		while (!localQueue.empty()){
//
//			Node_t node1 = localQueue.getNode();
////			zOpt = globalLB.load();
//			if (node1.ub < zOpt){
//				continue;
//				// pruned by bound.
//			}
//
//			auto result = explorer.process3(node1, zOpt);
//			nProcessed++;
//
//			if (result.success){
//				if (result.lb > zOpt){
//					// read lower bound.
//					zOpt = globalLB.load(memory_order_acquire);
//					if(result.lb > zOpt){
//						zOpt = result.lb;
//						globalLB.store(result.lb, memory_order_release); // check this in case of multiple threads.
//						const auto now = std::chrono::system_clock::now();
//						const auto t_c = std::chrono::system_clock::to_time_t(now);
//						cout << "thread: " << this_thread::get_id() << " , optimal LB: " << result.lb
//								<< " set at " << std::ctime(&t_c) << endl;
//					}
//				}
//				if (!result.nodes.empty()){
//					if (result.ub > zOpt){
//						localQueue.pushNodes(result.nodes);
//					}
//				}
//			}
//		}
//	}
//	cout << "Processed " << nProcessed << " nodes. " << endl;
//	explorer.displayCutStats();
//
//	//
//}

void DDSolver::startPThreadSolver() {
	// fill up work queue, make thread lock queue before accessing.
//	std::mutex queueLock;

	pair<CutContainer, CutContainer> cuts{CutContainer{FEASIBILITY}, CutContainer{OPTIMALITY}};

	NodeQueue tempQ1, tempQ2;
    vector<Node_t> t1, t2;
	vector<vector<Node_t>> vecvec{NUM_WORKERS};

	{

		// single iteration of node explorer.
		Node_t node {{},{}, numeric_limits<double>::lowest(), numeric_limits<double>::max(),0};
		NodeExplorer explorer {networkPtr};
		auto result = explorer.process3(node, numeric_limits<double>::lowest());

		if (result.lb > globalLB.load(memory_order_acquire))
			globalLB.store(result.lb, memory_order_release);
//		nodeQueue.pushNodes(result.nodes);
		int size = result.nodes.size();
		for (int i = 0; i < size; i++) {
		    if (i < size/2) {tempQ1.pushNode(result.nodes[i]); t1.push_back(result.nodes[i]);}
		    else {tempQ2.pushNode(result.nodes[i]); t2.push_back(result.nodes[i]);}
			vecvec[i%NUM_WORKERS].push_back(result.nodes[i]);
		}
		// vecvec[0] = t1;
		// vecvec[1] = t2;
		cuts = explorer.getCuts();

//		DDNode root{0};
//		root.globalLayer = 0;
//		DD restrictedDD{networkPtr, RESTRICTED};
//		auto cutset = restrictedDD.build(root);
//
//		// insert to cutset.
//		if (cutset)
//			nodeQueue.pushNodes(cutset.value());

		// cout << "Number of nodes in the queue: " << nodeQueue.size() << endl;

		// start thread
	}

    const unsigned nWorkers = 2;
	vector<thread> workerThreads{NUM_WORKERS};

    // vector<WorkerElement> workers(2);
	// workers[0].addNodesToWorker(vecvec[0]);
	// workers[1].addNodesToWorker(vecvec[1]);
	for (unsigned i = 0; i < NUM_WORKERS; i++) {
		workers[i].addNodesToWorker(vecvec[i]);
		workerThreads[i] = thread{&DDSolver::processWork3, this, i, cuts};
		// if (workerThreads[i].joinable()) workerThreads[i].join();
	}
    // workers[0].addNodes(t1); workers[1].addNodes(t2);
	std::thread master{&DDSolver::startMaster3, this};

	if (master.joinable()) master.join();
	for (unsigned i = 0; i < NUM_WORKERS; i++) {
		if (workerThreads[i].joinable()) workerThreads[i].join();
	}
	// std::thread worker{&DDSolver::processWork3, this, 0, cuts};
	// std::thread worker2{&DDSolver::processWork3, this, 1, cuts};
	// std::thread worker3{&DDSolver::processWork3, this, 2 , cuts};
	//
	// if (worker.joinable()) {
	// 	//cout << "Worker thread is joinable" << endl;
	// 	worker.join();
	// }
	// if (worker2.joinable()) worker2.join();
	// if (master.joinable()) master.join();
	// if (worker3.joinable()) worker3.join();
	//else {cout << "worker thread is not joinable" << endl;}

	//assert(nodeQueue.empty());
	//cout << "Number of nodes in the queue: " << nodeQueue.size() << endl;

	cout << "Optimal Solution : " << globalLB.load() << endl;

}


void DDSolver::processWork(unsigned int id, pair<CutContainer, CutContainer> cuts) {

	NodeExplorer explorer{networkPtr, cuts}; // get initial cuts later.
	auto& payload = workers[id];
	bool done = false;
	auto nodeVec  = payload.getNodes(done);
	NodeQueue localQueue{nodeVec}; // initialize localQueue later.

	double zOpt = globalLB.load(memory_order_acquire);
	size_t nProcessed  = 0;

	while (!isCompleted.load(memory_order_seq_cst)) {
		if (localQueue.empty()) {
			cout << "thread: " << this_thread::get_id() << " local queue is empty. indicating master." << endl;
			auto nodes = payload.getNodes(done);
			if (done) {cout << "thread: solver is finished" << endl; break;}
			localQueue.pushNodes(nodes);
		}
		while (!localQueue.empty()){
			Node_t node = localQueue.getNode();
			auto result = explorer.process3(node, zOpt);
			nProcessed++;
			if (result.success) {
				if (result.lb > zOpt){
					zOpt = globalLB.load(memory_order_acquire);
					if (result.lb > zOpt){
						globalLB.store(result.lb, memory_order_release);
						zOpt = result.lb;
						const auto now = std::chrono::system_clock::now();
						const auto t_c = std::chrono::system_clock::to_time_t(now);
						cout << "thread: " << this_thread::get_id() << " , optimal LB: " << result.lb
								<< " set at " << std::ctime(&t_c) << endl;
					}
				}

				if (result.ub > zOpt && !result.nodes.empty()){
					// add nodes to queue
					localQueue.pushNodes(result.nodes);
				}
			}

			//if master wants some work?
			scoped_lock l{payload.lock};
			auto st = payload.getStatus();
			if (st == MASTER_NEEDS_NODES) {
				// cout << "thread: " << this_thread::get_id() << " master needs nodes" << endl;
				auto n = localQueue.size();
				auto nodes = localQueue.getNodes(n/2);
				payload.addNodesToMaster(nodes); // payload status is updated here.
				cout << "thread " << this_thread::get_id() << " sent nodes to master" << endl;
			}
		}
	}
	// display stats.
	cout << "thread: " << this_thread::get_id() << " processed " << nProcessed << " nodes." << endl;
	explorer.displayCutStats();
}

void DDSolver::startMaster() {
	// master thread.
	const unsigned nWorkers = NUM_WORKERS;

	// constantly iterate
	while (true) {

		// look for idle workers and assign work.
		unsigned idleCount = 0;
		unsigned working = 0;

		for (unsigned i = 0; i < nWorkers; i++) {
			auto& worker = workers[i];
			auto status = worker.getStatus();

			if (status & STATUS::WORKER_NEEDS_NODES) {
				idleCount++;
			}
			else if (status & STATUS::WORKER_SHARED_NODES) {
				// take work from worker.
				auto nodes = worker.getNodesFromWorker();
				cout << "master received nodes from thread: " << i << endl;
				if (!nodes.empty())nodeQueue.pushNodes(nodes);
				working++;
				worker.setStatus(STATUS::MASTER_RECEIVED_NODES);
			}
			else if (status & STATUS::WORKER_WORKING) {
				// if
				working++;
			}
			else if (status & STATUS::NOT_ENOUGH_NODES_TO_SHARE) {
				// change the stauts

			}
		}

		if (idleCount == nWorkers &&  nodeQueue.empty()) {
			isCompleted.store(true, memory_order_seq_cst);
			// wake up all the workers.
			for (unsigned i = 0; i < nWorkers; i++) {
				auto& worker = workers[i];
				worker.setStatus(SOLVER_FINISHED);
				worker.cv.notify_one();
			}
			cout <<"Master: Solver is finished. Notified all workers" << endl;
			break;
		}

		if (idleCount > 0) {
			if (!nodeQueue.empty()) {
				// assign nodes to idle workers.
				int chunk = (nodeQueue.size()+idleCount)/idleCount;
				for (unsigned i = 0; i < nWorkers && !nodeQueue.empty() && idleCount; i++) {
					auto nodes = nodeQueue.getNodes(chunk);
					auto& worker = workers[i];
					// scoped_lock mut{worker.lock};
					if (worker.getStatus() & STATUS::WORKER_NEEDS_NODES) {
						worker.addNodesToWorker(nodes);
						idleCount--;
					}
					else {
						worker.askWorkerForNodes();
					}
				}

				if (idleCount > 0) {
					// ask nodes sequentially.
					for (unsigned i = 0; i < nWorkers; i++) {
						auto& worker = workers[i];
						if (worker.getStatus() & STATUS::WORKER_WORKING) {
							worker.askWorkerForNodes();

						}
					}

				}
			}
		}
	}

}
//
// uint8_t Payload::getStatus() const noexcept{
// 	return payloadStatus;
// 	// return status.load(memory_order_seq_cst);
// }
//
// void Payload::addNodes(vector<Node_t> nodes) { // called by master to insert nodes
// 	{
// 		// scoped_lock l{lock};
// 		nodes_ = nodes;
// 		payloadStatus = MASTER_ASSIGNED_NODES;
// 	}
// 	// status.store(MASTER_ASSIGNED_NODES, memory_order_seq_cst);
// 	cv.notify_one();
// }
//
// vector<Node_t> Payload::getNodes(bool& done) {
// 	// get nodes from master or payload.
//
// 	std::unique_lock l{lock};
// 	payloadStatus = WORKER_NEEDS_NODES;
// 	if (!nodes_.empty()) {
// 		/* either master placed some nodes initially, or worker placed some nodes previously
// 		   for master on master's request. Status flag should be zero. */
// 		auto nodes = move(nodes_);
// 		payloadStatus = WORKER_WORKING;
// 		// status.store(0, memory_order_release); // really necessary?
// 		return nodes;
// 	}
// 	// payloadStatus = WORKER_NEEDS_NODES;
// 	cv.wait(l, [&]{return ( payloadStatus == MASTER_ASSIGNED_NODES || payloadStatus == SOLVER_FINISHED);});
// 	if (payloadStatus == SOLVER_FINISHED) {cout <<"thread: finish status acknowledged"; done=true; return {};}
// 	payloadStatus = WORKER_WORKING;
// 	done = false;
// 	auto nodes = move(nodes_);
// 	string s = "Recevied " + to_string(nodes.size()) + " nodes from master.\n";
// 	cout << s;
// 	return nodes;
//
// 	// indicate master and wait.
// 	// payloadStatus = WORKER_NEEDS_NODES;
// 	// status.store(WORKER_NEEDS_NODES,memory_order_seq_cst);
//
//
// 	// while (true) {
// 	// 	auto st = status.load(memory_order_seq_cst);
// 	// 	std::unique_lock ul{lock};
// 	// 	cv.wait(ul, [&]{return !nodes_.empty() || (st == SOLVER_FINISHED);}); // add isCompleted flag check to break out.
// 	// 	if (nodes_.empty()) return{}; // master doesn't have any work to share. quit.
// 	// 	vector<Node_t> work = std::move(nodes_);
// 	// 	status.store(WORKER_WORKING, memory_order_seq_cst);
// 	// 	// nodes_.clear();
// 	// 	return work;
// 	// }
// }
//
// bool Payload::masterRequireNodes() const noexcept {
// 	// status.store(MASTER_NEEDS_NODES, memory_order_release);
// 	return  payloadStatus == MASTER_NEEDS_NODES;
// 	// return status.load(memory_order_relaxed) & MASTER_NEEDS_NODES;
// }
//
// void Payload::askWorkerForNodes() { // ask worker for more nodes.
// 	payloadStatus = MASTER_NEEDS_NODES;
// 	// status.store(MASTER_NEEDS_NODES, memory_order_seq_cst);
// }
//
//
// void Payload::addNodesToMaster(vector<Node_t> nodes) {
// 	{
// 		// scoped_lock l{lock};
// 		nodes_ = move(nodes);
// 		payloadStatus = WORKER_SHARED_NODES;
// 	}
// 	// status.store(WORKER_SHARED_NODES, memory_order_seq_cst);
// }
//
// vector<Node_t> Payload::getNodesFromWorker() {
// 	// return nodes from worker.
// 	vector<Node_t> nodes;
//
// 	{
// 		// scoped_lock l{lock};
// 		nodes = move(nodes_);
// 	}
// 	// set update flag
// 	// status.store(MASTER_RECEIVED_NODES, memory_order_seq_cst);
// 	payloadStatus = MASTER_RECEIVED_NODES;
// 	return nodes;
// }
//
// void Payload::setStatus(uint8_t status_) {
// 	// status.store(status_, memory_order_seq_cst);
// 	payloadStatus = status_;
// }

void DDSolver::startMaster2() {
	cout << "master starting" << endl;
	NodeQueue globalQueue;
	while (true) {
		//
		unsigned idle = 0, processing = 0;
		// iterate through every node and count idle workers
		for (unsigned i = 0; i < NUM_WORKERS; i++) {
			auto& worker = workers[i];
			scoped_lock l{worker.lock};

			// get status
			auto st = worker.getStatus();
			if (st == WORKER_WORKING) {
				processing++;
			}
			else if (st  == WORKER_NEEDS_NODES) idle++;
		}

		if (idle == NUM_WORKERS && globalQueue.empty()) {
			// solver is finished.
			isCompleted.store(true, memory_order_seq_cst);
			cout << "Master: Solver is finished. " << endl;
			for (int i = 0; i < NUM_WORKERS; i++) {
				auto& worker = workers[i];
				{
					scoped_lock l {worker.lock};
					worker.setStatus(SOLVER_FINISHED);
				}
				worker.cv.notify_one();
			}

			cout << "Master: Notified all workers" << endl;
			break; // exit from loop.
		}

		if (idle > 0){
			// ask for work.
			// cout << "Number of idle workers: " << idle << endl;
			for (unsigned i = 0; i < NUM_WORKERS; i++) {
				auto& worker = workers[i];
				scoped_lock l{worker.lock};
				auto st = worker.getStatus();
				if (st == STATUS::WORKER_WORKING) {
					worker.askWorkerForNodes(); // notify the working workers.
					string s = "Master: informed worker " + to_string(i) +" for nodes\n";
					cout << s;
				}
				else if ( st == STATUS::WORKER_NEEDS_NODES) {
					// share nodes
					if (!globalQueue.empty()) {
						// node queue is not empty.
						auto s = globalQueue.size();
						s /= 2;
						auto nodes = globalQueue.getNodes(s);
						worker.addNodesToWorker(nodes); // status update in the function.
						cout << "master: sent nodes to worker" << endl;
					}
				}
				else if (st == STATUS::WORKER_SHARED_NODES) {
					// TAKE nodes
					auto nodes = worker.getNodesFromWorker();
					string s = "Master: received " + to_string(nodes.size()) + " nodes\n";
					cout << s;
					globalQueue.pushNodes(nodes);
				}
			}
		}
	}
}

void DDSolver::processWork3(unsigned int id, pair<CutContainer, CutContainer> cuts) {

	NodeExplorer explorer{networkPtr, cuts};
	auto& payload = workers[id];
	bool done = false;
	auto nodeVec = payload.getNodes(done);
	NodeQueue localQueue{nodeVec};

	double zOpt = globalLB.load(memory_order_acquire);
	size_t nProcessed = 0;

	while (!isCompleted.load(memory_order_seq_cst)) {
		if (localQueue.empty()) {
			string s = "thread : " + to_string(id) + " local queue is empty. indicating master\n";
			cout << s;
			auto nodes = payload.getNodes(done);
			if (done) {cout<< "thread: solver is finished" << endl; break;}
			auto st = "thread " + to_string(id) + " received " + to_string(nodes.size()) + " nodes from master\n";
			cout << st;
			localQueue.pushNodes(nodes);
		}

		while (!localQueue.empty()) {
			Node_t node = localQueue.getNode();
			auto result = explorer.process3(node, zOpt);
			nProcessed++;
			if (result.success) {
				if (result.lb > zOpt) {
					zOpt = globalLB.load(memory_order_acquire);
					if (result.lb > zOpt) {
						globalLB.store(result.lb, memory_order_release);
						zOpt = result.lb;
						const auto now = std::chrono::system_clock::now();
						const auto t_c = std::chrono::system_clock::to_time_t(now);
						cout << "thread: " << id << " , optimal LB: " << result.lb
							<< " set at " << std::ctime(&t_c) << endl;
					}
				}

				if (result.ub > zOpt && !result.nodes.empty()) {
					localQueue.pushNodes(result.nodes);
				}
			}

			// if master wants nodes?
			{
				scoped_lock l{payload.lock};
				// auto st = payload.status;
				if (payload.status == MASTER_NEEDS_NODES) {
					auto n = localQueue.size();
					auto sz = static_cast<size_t>(ceil(n/2));
					auto nodes = localQueue.getNodes(sz);
					payload.nodes_ = move(nodes);
					payload.status = WORKER_SHARED_NODES;
					auto s = "thread " + to_string(id) + " sent " + to_string(payload.nodes_.size()) + " nodes  to master\n";
					cout<< s;
				}
			}
		}
	}

	auto s = "thread: " + to_string(id) + " processed " + to_string(nProcessed) + " nodes\n"; cout << s;
	explorer.displayCutStats();
}


vector<Node_t> Payload::getNodes(bool& done) {

	unique_lock l{lock};
	status = WORKER_NEEDS_NODES;
	if (!nodes_.empty()) {
		auto nodes = move(nodes_);
		status = WORKER_WORKING;
		return nodes;
	}
	cv.wait(l, [&]{return (status == MASTER_ASSIGNED_NODES|| status == SOLVER_FINISHED);});
	if (status == SOLVER_FINISHED) { done = true; return{};}
	// auto s = "Received " + to_string(nodes_.size()) + " nodes from master.\n"; cout << s;
	auto nodes = move(nodes_);
	done = false;
	status = WORKER_WORKING;
	return nodes;
}

void Payload::setStatus(uint8_t status_) { // assuming already locked.
	status = status_;
}

uint8_t Payload::getStatus() const noexcept {
	return status; // assuming already locked.
}

void Payload::addNodesToWorker(vector<Node_t> nodes) {
	// thread must be waiting for nodes.
	auto n = nodes.size();
	{
		scoped_lock l{lock};
		nodes_ = move(nodes);
		status = MASTER_ASSIGNED_NODES;
	}
	cv.notify_one(); // wake up worker.
	auto s = "Master sent " + to_string(n) + " nodes\n"; cout << s;
}

void Payload::askWorkerForNodes() {
	status = MASTER_NEEDS_NODES;
	// assuming already locked.
}

vector<Node_t> Payload::getNodesFromWorker() {
	// scoped_lock l{lock}; // assuming master locks before this function.
	if (status == WORKER_SHARED_NODES) {
		vector<Node_t> nodes;
		status = MASTER_RECEIVED_NODES;
		return move(nodes_);
	} // if worker has not enough nodes, ignore.
	return {};
}

void Payload::addNodesToMaster(vector<Node_t> nodes) {
	auto n = nodes.size(); // caller locked the mutex.
	nodes_ = move(nodes);
	status = WORKER_SHARED_NODES;
}


void DDSolver::startMaster3() {
	cout << "Starting master." << endl;
	NodeQueue globalQueue;

	while (true) {
		unsigned idle = 0, processing = 0;

		// iterate through all workers.
		for (unsigned i = 0; i < NUM_WORKERS; i++) {
			auto& worker = workers[i];
			bool added = false;
			{
				scoped_lock l{worker.lock};
				auto st = worker.status;
				if (st == WORKER_WORKING) {
					processing++;
				}
				else if (st == WORKER_SHARED_NODES) {
					// cout << "Some worker shared nodes" << endl;
					auto nodes = worker.nodes_;
					worker.nodes_.clear();
					worker.status = MASTER_RECEIVED_NODES;
					processing++;
					globalQueue.pushNodes(nodes);
					auto s= "Master received " + to_string(nodes.size()) +
						" nodes from worker: "+ to_string(i)+"\n" ; cout << s;
				}
				else if (st == WORKER_NEEDS_NODES) {
					if (!globalQueue.empty()) {
						// add nodes
						auto sz = globalQueue.size();
						sz  = static_cast<size_t> (ceil(static_cast<double>(sz)*0.5));
						auto nodes = globalQueue.getNodes(sz);
						worker.nodes_ = nodes; nodes.clear();
						worker.status = MASTER_ASSIGNED_NODES;
						added = true;
						auto s = "master added " + to_string(worker.nodes_.size()) +
							" nodes to worker " + to_string(i)+"\n"; cout <<s;
					}
					else idle++;
				}
				else if (st == NOT_ENOUGH_NODES_TO_SHARE || st == MASTER_RECEIVED_NODES) {
					// if not enough nodes to share or recevied nodes, change status.
					worker.status = WORKER_WORKING;
					processing++;
				}
				// else if (st == WORKER_SHARED_NODES) {
				// 	cout << "Some worker shared nodes" << endl;
				// 	// get nodes
				// 	auto nodes = worker.nodes_; worker.nodes_.clear();
				// 	worker.status = MASTER_RECEIVED_NODES;
				// 	processing++;
				// 	auto s = "Master received " + to_string(nodes.size()) +
				// 			" nodes from worker: "+ to_string(i)+"\n"; cout << s;
				// 	globalQueue.pushNodes(nodes);
				// }
			}
			if (added) worker.cv.notify_one();
		}

		if (idle == NUM_WORKERS && globalQueue.empty()) {
			// solver is finished.
			isCompleted.store(true, memory_order_seq_cst);
			cout << "Solver is finished" << endl;
			for (unsigned i = 0; i < NUM_WORKERS; i++) {
				auto& worker = workers[i];
				{
					scoped_lock l{worker.lock};
					worker.status = SOLVER_FINISHED;
				}
				worker.cv.notify_one();
			}
			cout << "master indicated all workers" <<endl;
			return;
		}

		if (idle > 0 ) {
			auto s = "Number of nodes in the queue: " + to_string(globalQueue.size()) + "\n";
			cout << s;
			// ask processing nodes for work.
			for (unsigned i = 0; i < NUM_WORKERS; i++) {
				auto& worker = workers[i];
				{
					scoped_lock l{worker.lock};
					if (worker.status == WORKER_WORKING) {
						// ask this worker.
						worker.status = MASTER_NEEDS_NODES;
					}
					else if (worker.status == NOT_ENOUGH_NODES_TO_SHARE || worker.status == MASTER_RECEIVED_NODES) {
						worker.status = WORKER_WORKING;
					}
				}
			}
		}
		this_thread::sleep_for(chrono::seconds(1));
	}
}
