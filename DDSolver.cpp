//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"
#include <random>
#include <chrono>
// #include <omp.h>
// #define OMP_NUM_THREADS 4
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

DDNode node2DDdfsNode(Node_t node) {
    DDNode newNode;
    newNode.states = set<int>(node.states.begin(), node.states.end());
    newNode.solutionVector = node.solutionVector;
    newNode.globalLayer = node.globalLayer;
    newNode.nodeLayer = 0;
    return newNode;
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
				if (payload.status == Payload::MASTER_NEEDS_NODES) {
					auto n = localQueue.size();
					auto sz = static_cast<size_t>(ceil(n/2));
					auto nodes = localQueue.getNodes(sz);
					payload.nodes_ = move(nodes);
					payload.status = Payload::WORKER_SHARED_NODES;
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
	{
		scoped_lock l{lock};
		nodes_ = move(nodes);
		// this request only occurs when the master specifically asked for nodes.
		// the previous payload status should be MASTER_NEEDS_NODES. any other status should not be appeared.
	}
	payloadStatus.store(WORKER_SHARED_NODES, memory_order::release);
	// nodes_ = move(nodes);
	// status = WORKER_SHARED_NODES;
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
				if (st == Payload::WORKER_WORKING) {
					processing++;
				}
				else if (st == Payload::WORKER_SHARED_NODES) {
					// cout << "Some worker shared nodes" << endl;
					auto nodes = worker.nodes_;
					worker.nodes_.clear();
					worker.status = Payload::MASTER_RECEIVED_NODES;
					processing++;
					globalQueue.pushNodes(nodes);
					auto s= "Master received " + to_string(nodes.size()) +
						" nodes from worker: "+ to_string(i)+"\n" ; cout << s;
				}
				else if (st == Payload::WORKER_NEEDS_NODES) {
					if (!globalQueue.empty()) {
						// add nodes
						auto sz = globalQueue.size();
						sz  = static_cast<size_t> (ceil(static_cast<double>(sz)*0.5));
						auto nodes = globalQueue.getNodes(sz);
						worker.nodes_ = nodes; nodes.clear();
						worker.status = Payload::MASTER_ASSIGNED_NODES;
						added = true;
						auto s = "master added " + to_string(worker.nodes_.size()) +
							" nodes to worker " + to_string(i)+"\n"; cout <<s;
					}
					else idle++;
				}
				else if (st == Payload::NOT_ENOUGH_NODES_TO_SHARE || st == Payload::MASTER_RECEIVED_NODES) {
					// if not enough nodes to share or recevied nodes, change status.
					worker.status = Payload::WORKER_WORKING;
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
					worker.status = Payload::SOLVER_FINISHED;
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
					if (worker.status == Payload::WORKER_WORKING) {
						// ask this worker.
						worker.status = Payload::MASTER_NEEDS_NODES;
					}
					else if (worker.status == Payload::NOT_ENOUGH_NODES_TO_SHARE || worker.status == Payload::MASTER_RECEIVED_NODES) {
						worker.status = Payload::WORKER_WORKING;
					}
				}
			}
		}
		this_thread::sleep_for(chrono::seconds(1));
	}
}


void DDSolver::Master::operator()(DDSolver &solver) {
	cout << "Master started" << endl;

	NodeQueue globalQueue;

	while (true) {
		unsigned idle = 0, processing = 0;

		for (unsigned i = 0; i < NUM_WORKERS; i++) {
			auto& worker = solver.workers[i];
			bool added = false;
			{
				scoped_lock l{worker.lock};
				auto st = worker.payloadStatus.load(memory_order::acquire);
				if (st == Payload::WORKER_WORKING) {
					processing++;
				}
				else if (st == Payload::WORKER_SHARED_NODES) {
					auto nodes = worker.nodes_;
					worker.nodes_.clear();
					worker.payloadStatus.store(Payload::MASTER_RECEIVED_NODES, memory_order::relaxed);
					processing++;
					globalQueue.pushNodes(nodes);
					auto s = "Master received " + to_string(worker.nodes_.size()) + " nodes from worker: " + to_string(i)+"\n" ; cout << s;
				}
				else if (st == Payload::WORKER_NEEDS_NODES) {
					if (!globalQueue.empty()) {
						// add ndoes
						auto sz = globalQueue.size();
						sz  = static_cast<size_t> (ceil(static_cast<double>(sz)*0.5));
						auto nodes = globalQueue.getNodes(sz);
						worker.nodes_ = move(nodes);
						added = true;
						worker.payloadStatus.store(Payload::STATUS::MASTER_ASSIGNED_NODES, memory_order::relaxed);
					}
					else idle++;
				}
				else if (st == Payload::STATUS::NOT_ENOUGH_NODES_TO_SHARE || st == Payload::STATUS::MASTER_RECEIVED_NODES) {
					worker.payloadStatus.store(Payload::STATUS::WORKER_WORKING, memory_order::relaxed);
					processing++;
				}

				if (worker.feasibilityCuts_ != nullptr) {
					fCutsGlobal.push_back(worker.feasibilityCuts_);
					worker.feasibilityCuts_ = nullptr;
				}
				if (worker.optimalityCuts_ != nullptr) {
					oCutsGlobal.push_back(worker.optimalityCuts_);
					worker.optimalityCuts_ = nullptr;
				}
			}
			if (added) worker.cv.notify_one();
		}
		if (idle == NUM_WORKERS && nodeQueue.empty()) {
			solver.isCompleted.store(true, memory_order::relaxed);
			for_each(solver.workers.begin(), solver.workers.end(), [](auto& worker) {
				worker.payloadStatus.store(Payload::STATUS::SOLVER_FINISHED, memory_order::release);
				worker.cv.notify_one();
			});
			cout << "Master indicated all the workers" << endl;
			return;
		}
		if (idle > 0) {
			// can I use memory_order_relax on all these updates and use single release order?
			for_each(solver.workers.begin(), solver.workers.end(), [](auto& worker) {
				worker.payloadStatus.compare_exchange_strong(Payload::STATUS::WORKER_WORKING, Payload::STATUS::MASTER_NEEDS_NODES, memory_order::acq_rel);
			});
		}
	}
}

void DDSolver::Worker::operator()(DDSolver &solver) {

	NodeExplorer explorer{solver.networkPtr};
	auto& payload = solver.workers[id];
	bool done = false;
	NodeQueue localQueue{payload.getNodes(done)};

	double zOpt = solver.globalLB.load(memory_order::relaxed);
	size_t nProcessed = 0;
	auto globalCuts = make_pair(fCutsGlobal, oCutsGlobal);

	uint prevCount = 0;

	while (!solver.isCompleted.load(memory_order::relaxed)) {
		if (localQueue.empty()) {
			string s = "Thread: " + to_string(id) + " local queue is empty.\n"; cout << s;
			auto nodes = payload.getNodes(done);
			if (done){ cout << ("Thread: "+to_string(id)+ " stopping...\n") << endl; break;}
			cout << ("Thread: "+to_string(id) + " received " + to_string(nodes.size()) +" nodes from master\n" ) <<endl;
			localQueue.pushNodes(nodes);
		}

		while (!localQueue.empty()) {
			Node_t node = localQueue.getNode();
			auto result = explorer.process4(node, zOpt,globalCuts);
			nProcessed++;
			// if cuts in the explorer exceeds limit? update to global.
			if ((explorer.optimalityCuts.cuts.size() + explorer.feasibilityCuts.cuts.size()) > LOCAL_CUTS_LIMIT) {
				shareCutsWithMaster(explorer, payload);
				// at this point explorer's cuts are empty.
			}
			if (solver.cutResources.getCount() > prevCount) { // TODO instead of polling continuously, poll periodically.
				// new cuts are added to the container, update local pointers.
				auto res = solver.cutResources.get(fCutsGlobal.size(), oCutsGlobal.size());
				if (!res.first.empty()) fCutsGlobal.insert(fCutsGlobal.end(), res.first.begin(), res.first.end());
				if (!res.second.empty()) oCutsGlobal.insert(oCutsGlobal.end(), res.second.begin(), res.second.end());
				prevCount = fCutsGlobal.size() + oCutsGlobal.size();
			}
			if (result.success) {
				if (is_poll_time(nProcessed) || result.lb > zOpt) { // poll for optimal value periodically.
					zOpt = solver.globalLB.load(memory_order::acquire);
					if (result.lb > zOpt) {
						while (!solver.globalLB.compare_exchange_strong(zOpt, result.lb, memory_order::release))
							{if(result.lb<zOpt) break;}
						if (zOpt == result.lb) {
							const auto now = std::chrono::system_clock::now();
							const auto t_c = std::chrono::system_clock::to_time_t(now);
							cout << "Thread: " << id << " , optimal LB: " << result.lb << " set at " << std::ctime(&t_c) << endl;
						}
						// share cuts with master.
					}
				}
				if (result.ub > zOpt && !result.nodes.empty()) localQueue.pushNodes(result.nodes);
			}
			// if master wants nodes?
			if (payload.payloadStatus.load(memory_order::relaxed) == Payload::STATUS::MASTER_NEEDS_NODES) {
				// scoped_lock l{payload.lock};
				auto n = localQueue.size();
				auto sz = static_cast<size_t>(ceil(static_cast<double>(n)/2));
				auto nodes = localQueue.getNodes(sz);
				payload.addNodesToMaster(nodes);
				shareCutsWithMaster(explorer, payload); // share cuts
			}
		}
	}
#ifdef SOLVER_STATS
	auto s = "Thread: " + to_string(id) + " processed " + to_string(nProcessed) + " nodes\n"; cout << s;
	explorer.displayCutStats();
#endif
}

void DDSolver::Worker::shareCutsWithMaster(NodeExplorer &explorer, Payload &payload) {
	// acquire lock
	// CutContainer* fcuts = new CutContainer(explorer.feasibilityCuts);
	// CutContainer* ocuts = new CutContainer(explorer.optimalityCuts);
	// clear clear cut containers.
	// explorer.feasibilityCuts.clearContainer();
	// explorer.optimalityCuts.clearContainer();
	{
		scoped_lock l{payload.lock};
		payload.feasibilityCuts_ = std::move(explorer.feasibilityCuts);
		payload.optimalityCuts_ = std::move(explorer.optimalityCuts);
		// atomic status update (with release order).
	}
}

void DDSolver::Master::addCutsToGlobal(DDSolver &solver) {

	// create cut containers defined in Inavap namespace.
	Inavap::CutContainer * fcuts = new Inavap::CutContainer();
	Inavap::CutContainer * ocuts = new Inavap::CutContainer();

	// insert fCuts to fcuts.
	for (auto container_p: fCutsGlobal) {
		for (auto cut : (*container_p).cuts) {
			Inavap::Cut c{cut.RHS, cut.cutCoeff};
			fcuts->insertCut(c);
		}
	}

	for (auto container_p : oCutsGlobal) {
		for (auto cut: container_p->cuts) {
			Inavap::Cut c{cut.RHS, cut.cutCoeff};
			ocuts->insertCut(c);
		}
	}
	// add to cutresource
	solver.cutResources.add({{fcuts}, {ocuts}});
}


