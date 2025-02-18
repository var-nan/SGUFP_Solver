//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"
#include <random>
#include <chrono>
// // #include <omp.h>
// // #define OMP_NUM_THREADS 4
// void DDSolver::NodeQueue::pushNodes(vector<Node_t> nodes) {
//     // push bunch of nodes.
//     // q.insert(q.end(), nodes.begin(), nodes.end());
//     for (auto& node : nodes) {
//         q.push(node);
//     }
// }
//
// void DDSolver::NodeQueue::pushNode(Node_t node) {
//     q.push(node);
// }
//
// Node_t DDSolver::NodeQueue::getNode() {
//     // auto node = q.back(); q.pop_back();
//     // auto node = q.front(); q.pop();
//     auto node = q.top(); q.pop();
//     return node;
// }
//
// vector<Node_t> DDSolver::NodeQueue::getNodes(size_t n = 8) {
//     vector<Node_t> nodes;
//
// 	while (!q.empty() && n) {
// 		//
// 		if (n--) {
// 			nodes.push_back(q.top()); q.pop();
// 		}
// 		// else break;
// 	}
//     //
//     // for (size_t i = 0; i < n; i++) {
//     //     nodes.push_back(q.back());
//     //     q.pop_back();
//     // }
//     return nodes;
// }
//
// Node_t DDSolver::getNode() {
//     return nodeQueue.getNode();
// }
//
// void DDSolver::initialize() {
//     // place root node to queue
//     Node_t node;
//
//     node.lb = std::numeric_limits<int64_t>::min();
//     node.ub = std::numeric_limits<int64_t>::max();
//     node.globalLayer = 0;
//
//     nodeQueue.pushNode(node);
// }
//
// void DDSolver::startSolve(optional<pair<CutContainer, CutContainer>> initialCuts = nullopt) {
//     if (initialCuts) {
//         NodeExplorer explorer{networkPtr, initialCuts.value()};
//         process(explorer);
//     }
//     else {
//         NodeExplorer explorer{networkPtr};
//         process(explorer);
//     }
// }
//
// DDNode node2DDdfsNode(Node_t node) {
//     DDNode newNode;
//     newNode.states = set<int>(node.states.begin(), node.states.end());
//     newNode.solutionVector = node.solutionVector;
//     newNode.globalLayer = node.globalLayer;
//     newNode.nodeLayer = 0;
//     return newNode;
// }
//
//
//
// void DDSolver::processWork3(unsigned int id, pair<CutContainer, CutContainer> cuts) {
//
// 	NodeExplorer explorer{networkPtr, cuts};
// 	auto& payload = workers[id];
// 	bool done = false;
// 	auto nodeVec = payload.getNodes(done);
// 	NodeQueue localQueue{nodeVec};
//
// 	double zOpt = globalLB.load(memory_order_acquire);
// 	size_t nProcessed = 0;
//
// 	while (!isCompleted.load(memory_order_seq_cst)) {
// 		if (localQueue.empty()) {
// 			string s = "thread : " + to_string(id) + " local queue is empty. indicating master\n";
// 			cout << s;
// 			auto nodes = payload.getNodes(done);
// 			if (done) {cout<< "thread: solver is finished" << endl; break;}
// 			auto st = "thread " + to_string(id) + " received " + to_string(nodes.size()) + " nodes from master\n";
// 			cout << st;
// 			localQueue.pushNodes(nodes);
// 		}
//
// 		while (!localQueue.empty()) {
// 			Node_t node = localQueue.getNode();
// 			auto result = explorer.process3(node, zOpt);
// 			nProcessed++;
// 			if (result.success) {
// 				if (result.lb > zOpt) {
// 					zOpt = globalLB.load(memory_order_acquire);
// 					if (result.lb > zOpt) {
// 						globalLB.store(result.lb, memory_order_release);
// 						zOpt = result.lb;
// 						const auto now = std::chrono::system_clock::now();
// 						const auto t_c = std::chrono::system_clock::to_time_t(now);
// 						cout << "thread: " << id << " , optimal LB: " << result.lb
// 							<< " set at " << std::ctime(&t_c) << endl;
// 					}
// 				}
//
// 				if (result.ub > zOpt && !result.nodes.empty()) {
// 					localQueue.pushNodes(result.nodes);
// 				}
// 			}
//
// 			// if master wants nodes?
// 			{
// 				scoped_lock l{payload.lock};
// 				// auto st = payload.status;
// 				if (payload.status == Payload::MASTER_NEEDS_NODES) {
// 					auto n = localQueue.size();
// 					auto sz = static_cast<size_t>(ceil(n/2));
// 					auto nodes = localQueue.getNodes(sz);
// 					payload.nodes_ = move(nodes);
// 					payload.status = Payload::WORKER_SHARED_NODES;
// 					auto s = "thread " + to_string(id) + " sent " + to_string(payload.nodes_.size()) + " nodes  to master\n";
// 					cout<< s;
// 				}
// 			}
// 		}
// 	}
//
// 	auto s = "thread: " + to_string(id) + " processed " + to_string(nProcessed) + " nodes\n"; cout << s;
// 	explorer.displayCutStats();
// }
//
//
// vector<Node_t> Payload::getNodes(bool& done) {
//
// 	unique_lock l{lock};
// 	status = WORKER_NEEDS_NODES;
// 	if (!nodes_.empty()) {
// 		auto nodes = move(nodes_);
// 		status = WORKER_WORKING;
// 		return nodes;
// 	}
// 	cv.wait(l, [&]{return (status == MASTER_ASSIGNED_NODES|| status == SOLVER_FINISHED);});
// 	if (status == SOLVER_FINISHED) { done = true; return{};}
// 	// auto s = "Received " + to_string(nodes_.size()) + " nodes from master.\n"; cout << s;
// 	auto nodes = move(nodes_);
// 	done = false;
// 	status = WORKER_WORKING;
// 	return nodes;
// }
//
// void Payload::setStatus(uint8_t status_) { // assuming already locked.
// 	status = status_;
// }
//
// uint8_t Payload::getStatus() const noexcept {
// 	return status; // assuming already locked.
// }
//
// void Payload::addNodesToWorker(vector<Node_t> nodes) {
// 	// thread must be waiting for nodes.
// 	auto n = nodes.size();
// 	{
// 		scoped_lock l{lock};
// 		nodes_ = move(nodes);
// 		status = MASTER_ASSIGNED_NODES;
// 	}
// 	cv.notify_one(); // wake up worker.
// 	auto s = "Master sent " + to_string(n) + " nodes\n"; cout << s;
// }
//
// void Payload::askWorkerForNodes() {
// 	status = MASTER_NEEDS_NODES;
// 	// assuming already locked.
// }
//
// vector<Node_t> Payload::getNodesFromWorker() {
// 	// scoped_lock l{lock}; // assuming master locks before this function.
// 	if (status == WORKER_SHARED_NODES) {
// 		vector<Node_t> nodes;
// 		status = MASTER_RECEIVED_NODES;
// 		return move(nodes_);
// 	} // if worker has not enough nodes, ignore.
// 	return {};
// }
//
// void Payload::addNodesToMaster(vector<Node_t> nodes) {
// 	auto n = nodes.size(); // caller locked the mutex.
// 	{
// 		scoped_lock l{lock};
// 		nodes_ = move(nodes);
// 		// this request only occurs when the master specifically asked for nodes.
// 		// the previous payload status should be MASTER_NEEDS_NODES. any other status should not be appeared.
// 	}
// 	payloadStatus.store(WORKER_SHARED_NODES, memory_order::release);
// 	// nodes_ = move(nodes);
// 	// status = WORKER_SHARED_NODES;
// }
//
//
// void DDSolver::startMaster3() {
// 	cout << "Starting master." << endl;
// 	NodeQueue globalQueue;
//
// 	while (true) {
// 		unsigned idle = 0, processing = 0;
//
// 		// iterate through all workers.
// 		for (unsigned i = 0; i < NUM_WORKERS; i++) {
// 			auto& worker = workers[i];
// 			bool added = false;
// 			{
// 				scoped_lock l{worker.lock};
// 				auto st = worker.status;
// 				if (st == Payload::WORKER_WORKING) {
// 					processing++;
// 				}
// 				else if (st == Payload::WORKER_SHARED_NODES) {
// 					// cout << "Some worker shared nodes" << endl;
// 					auto nodes = worker.nodes_;
// 					worker.nodes_.clear();
// 					worker.status = Payload::MASTER_RECEIVED_NODES;
// 					processing++;
// 					globalQueue.pushNodes(nodes);
// 					auto s= "Master received " + to_string(nodes.size()) +
// 						" nodes from worker: "+ to_string(i)+"\n" ; cout << s;
// 				}
// 				else if (st == Payload::WORKER_NEEDS_NODES) {
// 					if (!globalQueue.empty()) {
// 						// add nodes
// 						auto sz = globalQueue.size();
// 						sz  = static_cast<size_t> (ceil(static_cast<double>(sz)*0.5));
// 						auto nodes = globalQueue.getNodes(sz);
// 						worker.nodes_ = nodes; nodes.clear();
// 						worker.status = Payload::MASTER_ASSIGNED_NODES;
// 						added = true;
// 						auto s = "master added " + to_string(worker.nodes_.size()) +
// 							" nodes to worker " + to_string(i)+"\n"; cout <<s;
// 					}
// 					else idle++;
// 				}
// 				else if (st == Payload::NOT_ENOUGH_NODES_TO_SHARE || st == Payload::MASTER_RECEIVED_NODES) {
// 					// if not enough nodes to share or recevied nodes, change status.
// 					worker.status = Payload::WORKER_WORKING;
// 					processing++;
// 				}
// 				// else if (st == WORKER_SHARED_NODES) {
// 				// 	cout << "Some worker shared nodes" << endl;
// 				// 	// get nodes
// 				// 	auto nodes = worker.nodes_; worker.nodes_.clear();
// 				// 	worker.status = MASTER_RECEIVED_NODES;
// 				// 	processing++;
// 				// 	auto s = "Master received " + to_string(nodes.size()) +
// 				// 			" nodes from worker: "+ to_string(i)+"\n"; cout << s;
// 				// 	globalQueue.pushNodes(nodes);
// 				// }
// 			}
// 			if (added) worker.cv.notify_one();
// 		}
//
// 		if (idle == NUM_WORKERS && globalQueue.empty()) {
// 			// solver is finished.
// 			isCompleted.store(true, memory_order_seq_cst);
// 			cout << "Solver is finished" << endl;
// 			for (unsigned i = 0; i < NUM_WORKERS; i++) {
// 				auto& worker = workers[i];
// 				{
// 					scoped_lock l{worker.lock};
// 					worker.status = Payload::SOLVER_FINISHED;
// 				}
// 				worker.cv.notify_one();
// 			}
// 			cout << "master indicated all workers" <<endl;
// 			return;
// 		}
//
// 		if (idle > 0 ) {
// 			auto s = "Number of nodes in the queue: " + to_string(globalQueue.size()) + "\n";
// 			cout << s;
// 			// ask processing nodes for work.
// 			for (unsigned i = 0; i < NUM_WORKERS; i++) {
// 				auto& worker = workers[i];
// 				{
// 					scoped_lock l{worker.lock};
// 					if (worker.status == Payload::WORKER_WORKING) {
// 						// ask this worker.
// 						worker.status = Payload::MASTER_NEEDS_NODES;
// 					}
// 					else if (worker.status == Payload::NOT_ENOUGH_NODES_TO_SHARE || worker.status == Payload::MASTER_RECEIVED_NODES) {
// 						worker.status = Payload::WORKER_WORKING;
// 					}
// 				}
// 			}
// 		}
// 		this_thread::sleep_for(chrono::seconds(1));
// 	}
// }

//
// void DDSolver::Master::operator()(DDSolver &solver) {
// 	cout << "Master started" << endl;
//
// 	NodeQueue globalQueue;
//
// 	while (true) {
// 		unsigned idle = 0, processing = 0;
//
// 		for (unsigned i = 0; i < NUM_WORKERS; i++) {
// 			auto& worker = solver.workers[i];
// 			bool added = false;
// 			{
// 				scoped_lock l{worker.lock};
// 				auto st = worker.payloadStatus.load(memory_order::acquire);
// 				if (st == Payload::WORKER_WORKING) {
// 					processing++;
// 				}
// 				else if (st == Payload::WORKER_SHARED_NODES) {
// 					auto nodes = worker.nodes_;
// 					worker.nodes_.clear();
// 					worker.payloadStatus.store(Payload::MASTER_RECEIVED_NODES, memory_order::relaxed);
// 					processing++;
// 					globalQueue.pushNodes(nodes);
// 					auto s = "Master received " + to_string(worker.nodes_.size()) + " nodes from worker: " + to_string(i)+"\n" ; cout << s;
// 				}
// 				else if (st == Payload::WORKER_NEEDS_NODES) {
// 					if (!globalQueue.empty()) {
// 						// add ndoes
// 						auto sz = globalQueue.size();
// 						sz  = static_cast<size_t> (ceil(static_cast<double>(sz)*0.5));
// 						auto nodes = globalQueue.getNodes(sz);
// 						worker.nodes_ = move(nodes);
// 						added = true;
// 						worker.payloadStatus.store(Payload::STATUS::MASTER_ASSIGNED_NODES, memory_order::relaxed);
// 					}
// 					else idle++;
// 				}
// 				else if (st == Payload::STATUS::NOT_ENOUGH_NODES_TO_SHARE || st == Payload::STATUS::MASTER_RECEIVED_NODES) {
// 					worker.payloadStatus.store(Payload::STATUS::WORKER_WORKING, memory_order::relaxed);
// 					processing++;
// 				}
//
// 				if (worker.feasibilityCuts_ != nullptr) {
// 					fCutsGlobal.push_back(worker.feasibilityCuts_);
// 					worker.feasibilityCuts_ = nullptr;
// 				}
// 				if (worker.optimalityCuts_ != nullptr) {
// 					oCutsGlobal.push_back(worker.optimalityCuts_);
// 					worker.optimalityCuts_ = nullptr;
// 				}
// 			}
// 			if (added) worker.cv.notify_one();
// 		}
// 		if (idle == NUM_WORKERS && nodeQueue.empty()) {
// 			solver.isCompleted.store(true, memory_order::relaxed);
// 			for_each(solver.workers.begin(), solver.workers.end(), [](auto& worker) {
// 				worker.payloadStatus.store(Payload::STATUS::SOLVER_FINISHED, memory_order::release);
// 				worker.cv.notify_one();
// 			});
// 			cout << "Master indicated all the workers" << endl;
// 			return;
// 		}
// 		if (idle > 0) {
// 			// can I use memory_order_relax on all these updates and use single release order?
// 			for_each(solver.workers.begin(), solver.workers.end(), [](auto& worker) {
// 				worker.payloadStatus.compare_exchange_strong(Payload::STATUS::WORKER_WORKING, Payload::STATUS::MASTER_NEEDS_NODES, memory_order::acq_rel);
// 			});
// 		}
// 	}
// }
//
// void DDSolver::Worker::operator()(DDSolver &solver) {
//
// 	NodeExplorer explorer{solver.networkPtr};
// 	auto& payload = solver.workers[id];
// 	bool done = false;
// 	NodeQueue localQueue{payload.getNodes(done)};
//
// 	double zOpt = solver.globalLB.load(memory_order::relaxed);
// 	size_t nProcessed = 0;
// 	auto globalCuts = make_pair(fCutsGlobal, oCutsGlobal);
//
// 	uint prevCount = 0;
//
// 	while (!solver.isCompleted.load(memory_order::relaxed)) {
// 		if (localQueue.empty()) {
// 			string s = "Thread: " + to_string(id) + " local queue is empty.\n"; cout << s;
// 			auto nodes = payload.getNodes(done);
// 			if (done){ cout << ("Thread: "+to_string(id)+ " stopping...\n") << endl; break;}
// 			cout << ("Thread: "+to_string(id) + " received " + to_string(nodes.size()) +" nodes from master\n" ) <<endl;
// 			localQueue.pushNodes(nodes);
// 		}
//
// 		while (!localQueue.empty()) {
// 			Node_t node = localQueue.getNode();
// 			auto result = explorer.process4(node, zOpt,globalCuts);
// 			nProcessed++;
// 			// if cuts in the explorer exceeds limit? update to global.
// 			if ((explorer.optimalityCuts.cuts.size() + explorer.feasibilityCuts.cuts.size()) > LOCAL_CUTS_LIMIT) {
// 				shareCutsWithMaster(explorer, payload);
// 				// at this point explorer's cuts are empty.
// 			}
// 			if (solver.cutResources.getCount() > prevCount) { // TODO instead of polling continuously, poll periodically.
// 				// new cuts are added to the container, update local pointers.
// 				auto res = solver.cutResources.get(fCutsGlobal.size(), oCutsGlobal.size());
// 				if (!res.first.empty()) fCutsGlobal.insert(fCutsGlobal.end(), res.first.begin(), res.first.end());
// 				if (!res.second.empty()) oCutsGlobal.insert(oCutsGlobal.end(), res.second.begin(), res.second.end());
// 				prevCount = fCutsGlobal.size() + oCutsGlobal.size();
// 			}
// 			if (result.success) {
// 				if (is_poll_time(nProcessed) || result.lb > zOpt) { // poll for optimal value periodically.
// 					zOpt = solver.globalLB.load(memory_order::acquire);
// 					if (result.lb > zOpt) {
// 						while (!solver.globalLB.compare_exchange_strong(zOpt, result.lb, memory_order::release))
// 							{if(result.lb<zOpt) break;}
// 						if (zOpt == result.lb) {
// 							const auto now = std::chrono::system_clock::now();
// 							const auto t_c = std::chrono::system_clock::to_time_t(now);
// 							cout << "Thread: " << id << " , optimal LB: " << result.lb << " set at " << std::ctime(&t_c) << endl;
// 						}
// 						// share cuts with master.
// 					}
// 				}
// 				if (result.ub > zOpt && !result.nodes.empty()) localQueue.pushNodes(result.nodes);
// 			}
// 			// if master wants nodes?
// 			if (payload.payloadStatus.load(memory_order::relaxed) == Payload::STATUS::MASTER_NEEDS_NODES) {
// 				// scoped_lock l{payload.lock};
// 				auto n = localQueue.size();
// 				auto sz = static_cast<size_t>(ceil(static_cast<double>(n)/2));
// 				auto nodes = localQueue.getNodes(sz);
// 				payload.addNodesToMaster(nodes);
// 				shareCutsWithMaster(explorer, payload); // share cuts
// 			}
// 		}
// 	}
// #ifdef SOLVER_STATS
// 	auto s = "Thread: " + to_string(id) + " processed " + to_string(nProcessed) + " nodes\n"; cout << s;
// 	explorer.displayCutStats();
// #endif
// }
//
// void DDSolver::Worker::shareCutsWithMaster(NodeExplorer &explorer, Payload &payload) {
// 	// acquire lock
// 	// CutContainer* fcuts = new CutContainer(explorer.feasibilityCuts);
// 	// CutContainer* ocuts = new CutContainer(explorer.optimalityCuts);
// 	// clear clear cut containers.
// 	// explorer.feasibilityCuts.clearContainer();
// 	// explorer.optimalityCuts.clearContainer();
// 	{
// 		scoped_lock l{payload.lock};
// 		payload.feasibilityCuts_ = std::move(explorer.feasibilityCuts);
// 		payload.optimalityCuts_ = std::move(explorer.optimalityCuts);
// 		// atomic status update (with release order).
// 	}
// }
//
// void DDSolver::Master::addCutsToGlobal(DDSolver &solver) {
//
// 	// create cut containers defined in Inavap namespace.
// 	Inavap::CutContainer * fcuts = new Inavap::CutContainer();
// 	Inavap::CutContainer * ocuts = new Inavap::CutContainer();
//
// 	// insert fCuts to fcuts.
// 	for (auto container_p: fCutsGlobal) {
// 		for (auto cut : (*container_p).cuts) {
// 			Inavap::Cut c{cut.RHS, cut.cutCoeff};
// 			fcuts->insertCut(c);
// 		}
// 	}
//
// 	for (auto container_p : oCutsGlobal) {
// 		for (auto cut: container_p->cuts) {
// 			Inavap::Cut c{cut.RHS, cut.cutCoeff};
// 			ocuts->insertCut(c);
// 		}
// 	}
// 	// add to cutresource
// 	solver.cutResources.add({{fcuts}, {ocuts}});
// }


// bool Inavap::DDSolver::Master::isCutsShared(Payload &worker) {
// 	// check if worker shared cuts?
// }

/**
 * Adds the given Cuts to the global space.
 * @param solver
 * @param cuts
 */
void Inavap::DDSolver::Master::addCutsToGlobal(DDSolver &solver, pair<vector<Inavap::CutContainer>, vector<Inavap::CutContainer>> &cuts) {

	vector<CutContainer *> feasibilityPointers, optimalityPointers;

	size_t nfCuts = accumulate(cuts.first.begin(), cuts.first.end(), 0,
		[](size_t sum, const CutContainer &cont) { return sum + cont.size(); });
	// size_t nfCuts = 0;
	// for (const auto &cutContainer : cuts.first) {
	// 	nfCuts += cutContainer.size();
	// }

	// TODO: add master's local cuts to the global.

	if (nfCuts > FEASIBILITY_CONTAINER_CAPACITY) {
		uint nContainers = static_cast<uint> (floor(nfCuts/FEASIBILITY_CONTAINER_CAPACITY));
		if (nContainers) {
			feasibilityPointers = vector<CutContainer *>(nContainers);
			uint c_index = 0;
			uint n_so_far = 0;
			for (auto& container : cuts.first) {
				for (auto&& cut : container) { // don't know if it causes problem?
					feasibilityPointers[c_index]->insertCut(std::move(cut));
					n_so_far++;
					if (n_so_far == FEASIBILITY_CONTAINER_CAPACITY) c_index++;
				}
			}
		}
	}

	size_t nOCuts = accumulate(cuts.second.begin(), cuts.second.end(), 0,
		[](size_t sum, const CutContainer &cont) { return sum + cont.size(); });
	if (nOCuts > OPTIMALITY_CONTAINER_CAPACITY) {
		uint nContainers = static_cast<uint>(floor(nOCuts/OPTIMALITY_CONTAINER_CAPACITY));
		if (nContainers) {
			optimalityPointers = vector<CutContainer *> (nContainers);
			uint c_index = 0;
			uint n_so_far = 0;
			for (auto& container : cuts.second) {
				for (auto&& cut: container) { // capturing by references (without copying).
					optimalityPointers[c_index]->insertCut(std::move(cut));
					n_so_far++;
					if (n_so_far == OPTIMALITY_CONTAINER_CAPACITY) c_index++;
				}
			}
		}
	}

	solver.CutResources.add(make_pair(feasibilityPointers, optimalityPointers));
	// if (nOCuts > OPTIMALITY_CONTAINER_CAPACITY) CutContainer
}

/**
 * Master thread enters to worker mode and process nodes from its local queue.
 * @param solver
 * @param explorer - Node explorer of master thread.
 * @param n - Number of nodes to process.
 * @return - number of nodes processed by the master in this function call.
 */
size_t Inavap::DDSolver::Master::processNodes(Inavap::DDSolver &solver, Inavap::NodeExplorer &explorer, size_t n) {
	// worker-mode.
	double zOpt = solver.optimal;
	size_t nProcessed = 0;

	while (!nodeQueue.empty() && n--) {
		auto result = explorer.process(nodeQueue.popNode(), zOpt, feasCutsGlobal, optCutsGlobal);
		nProcessed++;
		// TODO: what to do with local cuts?
		if (result.status) {
			if (result.lb > zOpt) {
				// use separate lock?
				omp_set_lock(&solver.optimal_lock);
				if (result.lb > solver.optimal)
					solver.optimal = result.lb;
				zOpt = solver.optimal;
				omp_unset_lock(&solver.optimal_lock);
			}
			if (result.ub > zOpt && !result.nodes.empty()) nodeQueue.pushNodes(result.nodes);
		}
	}
	return nProcessed;
}

void Inavap::DDSolver::Master::startMaster(DDSolver &solver) {
	cout << "Starting Master thread" << endl;

	NodeExplorer explorer{networkPtr};
	size_t nProcessed = 0;

	// check if the worker shared any cuts?
	auto func = [&](Payload &worker) { // pair of feasibility, optimality cuts.
		// need to lock again?
		pair<CutContainer, CutContainer> cuts;
		if (!worker.fCuts.empty()) {
			cuts.first = std::move(worker.fCuts);
		}
		if (!worker.oCuts.empty()) {
			cuts.second = std::move(worker.oCuts);
		}
		return cuts;
	};

	while (!solver.isCompleted) {
		// number of idle workers, and number of processing workers.
		uint16_t idle = 0, processing = 0;

		// containers to hold the cuts shared by the workers.
		vector<CutContainer> feasibilityCuts;
		vector<CutContainer> optimalityCuts;

		// iterate through all workers.
		for (uint16_t i = 0; i < solver.N_WORKERS; i++) {
			auto& worker = solver.payloads[i];
			pair<CutContainer, CutContainer> workerCuts;
			omp_set_lock(&worker.lock);
			{
				/* Worker might share cuts anytime. Instead of double locking for cuts and status check, acquire lock
				 * on entire payload before operation starts.*/
				// workerCuts = func(worker);
				if (worker.payloadStatus == Payload::WORKER_WORKING) {
					processing++;
				}
				else if (worker.payloadStatus == Payload::WORKER_SHARED_NODES) {
					auto message = "Thread " +std::to_string(worker.id) + " sent " + std::to_string(worker.payloadNodes.size()) + " nodes."; cout << message << endl;
                    auto nodes = std::move(worker.payloadNodes);
                    worker.payloadStatus = Payload::MASTER_RECEIVED_NODES;
                    processing++;
                    nodeQueue.pushNodes(nodes);
				}
				else if (worker.payloadStatus == Payload::WORKER_NEEDS_NODES) {
					// status cannot be changed by the worker once reached here.
					if(!nodeQueue.empty()) {
						{ //
							uint size = nodeQueue.size();
							size = static_cast<size_t> (ceil(size*PROPORTION_OF_SHARE));
							auto nodes = nodeQueue.popNodes(size);
							worker.payloadNodes = std::move(nodes);
							worker.payloadStatus = Payload::MASTER_ASSIGNED_NODES;
							auto msg = "Master sent " + to_string(worker.payloadNodes.size()) + " nodes to worker " + std::to_string(worker.id); cout << msg << endl;
						}
						processing++;
						// worker.cv.notify_one(); // signal after lock release instead.
						// watchStatus = 1;
					}
					else idle++;
				}
				else if (worker.payloadStatus == Payload::NOT_ENOUGH_NODES_TO_SHARE ||
							worker.payloadStatus == Payload::MASTER_RECEIVED_NODES) {
					worker.payloadStatus = Payload::WORKER_WORKING;
					processing++;
				}
			}
			omp_unset_lock(&worker.lock);
			// if (watchStatus) worker.cv.notify_one(); // since no conditional variable, remove this.
			// add cuts to global.
			if (!workerCuts.first.empty()) feasibilityCuts.push_back( std::move(workerCuts.first));
			if (!workerCuts.second.empty()) optimalityCuts.push_back(std::move(workerCuts.second));
		}

		if (idle == solver.N_WORKERS && nodeQueue.empty()) {
			cout << "Indicating all workers" << endl;
			omp_set_lock(&solver.bool_lock);
			solver.isCompleted = true;
			omp_unset_lock(&solver.bool_lock);
			// all workers must be in spin  loop, eventually reads the status flag from payload.
			for (auto& worker : solver.payloads) {
				omp_set_lock(&worker.lock);
				worker.payloadStatus = Payload::SOLVER_FINISHED;
				omp_unset_lock(&worker.lock);
			}
			cout << "Solver finished. "<< endl;
			return;
		}

		if (idle > 0 && false) { // do not ask nodes from workers.
			// some workers are idle, get nodes from busy workers.
			for_each(solver.payloads.begin(), solver.payloads.end(), [](auto &worker) {
				omp_set_lock(&worker.lock);
				if (worker.payloadStatus == Payload::WORKER_WORKING)
					worker.payloadStatus = Payload::MASTER_NEEDS_NODES;
				else if (worker.payloadStatus == Payload::NOT_ENOUGH_NODES_TO_SHARE ||
						worker.payloadStatus == Payload::MASTER_RECEIVED_NODES)
					worker.payloadStatus = Payload::WORKER_WORKING;
					// either this worker don't have enough nodes to share or it shared earlier.
			});
		}

		// thread sleep?
		// update local cuts to global.
		auto p = make_pair(feasibilityCuts, optimalityCuts);
		addCutsToGlobal(solver, p);

		sleep(2);
		// process some nodes from local queue?
		// nProcessed += processNodes(solver, explorer, 4);
	}
	// post-completion tasks by master?
}

void Inavap::DDSolver::Worker::startWorker(DDSolver *solver) {
	// change to pointer of DDSolver.

	/* function: returns true if it is time to read the atomic variable of global optimal, else increases counter. */
	auto is_poll_time = [](size_t &counter) { // resets counter when reaches limit.
       if (counter == POLL_FREQUENCY) { counter = 0; return true; }
        counter++; return false;
    };

	/* function: returns cuts from the node explorer */
    // auto f_shareCuts = [](NodeExplorer &explorer) { // share cuts with master if exceeds limit.
    //     pair<optional<CutContainer>, optional<CutContainer>> cuts;
    //     if (explorer.feasibilityCuts.size() > FEASIBILITY_CONTAINER_LIMIT )
    //         cuts.first.value() = move(explorer.feasibilityCuts);
    //     if (explorer.optimalityCuts.size() > OPTIMALITY_CONTAINER_LIMIT)
    //         cuts.second.value() = move(explorer.optimalityCuts);
    //     return cuts;
    // };


	NodeExplorer explorer{networkPtr};
	auto& payload = solver->payloads[id];
	payload.id = id;
	// auto str = "Worker: " + to_string(id) + " , nodes: " + to_string(payload.nodes.size()) + "\n"; cout << str <<  endl;
	NodeQueue localQueue;
	double zOpt = solver->optimal;

	uint8_t done = 0; // flag to indicate solver is finished.
	size_t nProcessed = 0; // # nodes processed so far.
	size_t counter = 0;
	uint prevCount = 0; // # newly generated cuts since last shared with global.


	while (!solver->isSolverFinsished()) {

		if (localQueue.empty()) {
			// indicate master for nodes. wait until master responds.
			auto nodes = payload.getNodes(done);
			if (done) { cout << "finished flag received"  << endl;break;} // solver is finished. isFinished flag has been set?
			localQueue.pushNodes(nodes);
		}

		while (!localQueue.empty()) {
			if (omp_get_thread_num() == 1) {
				if (zOpt >= -2190) {
					cout << "# nodes in queue : " << localQueue.size() << endl;
				}
			}
			// cout << "worker: " << id << " processing " << endl;

			// get updated cuts from global space.
			// if (solver->CutResources.getCount() > prevCount) { // new cuts are added to global space.
			// 	auto res = solver->CutResources.get(feasCutsGlobal.size(), optCutsGlobal.size());
			// 	if (!res.first.empty()) feasCutsGlobal.insert(feasCutsGlobal.end(), res.first.begin(), res.first.end());
			// 	if (!res.second.empty()) optCutsGlobal.insert(optCutsGlobal.end(), res.second.begin(), res.second.end());
			// 	prevCount = feasCutsGlobal.size() + optCutsGlobal.size();
			// }

			Node node = localQueue.popNode();
			auto result = explorer.process(node, zOpt, feasCutsGlobal, optCutsGlobal);
			// auto result = explorer.process2(localQueue.popNode(), zOpt);
			nProcessed++;
			// auto str = "Worker : " + to_string(id) + " processed a node\n"; cout << str << endl;

			if (result.status) {
				//
				if (result.lb > zOpt) {
					// porting to OpenMP.
					omp_set_lock(&solver->optimal_lock);
					if (result.lb > solver->optimal)
						solver->optimal = result.lb;
					zOpt = solver->optimal;
					omp_unset_lock(&solver->optimal_lock);

					 if (result.lb == zOpt) {
					 	const auto now = std::chrono::system_clock::now();
					 	const auto t_c = std::chrono::system_clock::to_time_t(now);
					 	cout << "Thread: " << id << " , optimal lb: " << zOpt << " set at, " << std::ctime(&t_c) << endl;
					 }

					counter = 0;
				}
				else if (is_poll_time(counter)) { // check periodically for latest optimal value.
					omp_set_lock(&solver->optimal_lock);
					zOpt = solver->optimal;
					omp_unset_lock(&solver->optimal_lock);
				}
				if (result.ub > zOpt && !result.nodes.empty())
					localQueue.pushNodes(result.nodes);
			}

			// look into explorer's local cuts.
			// if ((explorer.feasibilityCuts.size() + explorer.optimalityCuts.size()) > LOCAL_CUTS_LIMIT) {
			// 	// use different strategy.
			// 	scoped_lock l{payload.lock};
			// 	if (explorer.feasibilityCuts.size() > FEASIBILITY_CONTAINER_LIMIT) payload.fCuts = std::move(explorer.feasibilityCuts);
			// 	if (explorer.optimalityCuts.size() > OPTIMALITY_CONTAINER_LIMIT) payload.oCuts = std::move(explorer.optimalityCuts);
			// 	// shareCutsWithMaster(explorer, payload);
			// 	// At this point, node explorer cuts are empty.
			// }

            // if payload status is 'MASTER NEEDS NODES', then master cannot change it.
			// TODO: instead of checking in every iteration, check periodically
			// porting to OpenMP
			// omp_set_lock(&payload.lock);
   //          if (payload.payloadStatus == Payload::MASTER_NEEDS_NODES) {
   //              uint n = localQueue.size();
   //              auto sz = static_cast<size_t> (ceil(n*0.4));
   //              auto nodes = localQueue.popNodes(sz);
   //          	/* Do not share cuts with master. */
   //          	// auto [fst, snd] = f_shareCuts(explorer);
   //              // scoped_lock l{payload.lock};
   //              payload.nodes = std::move(nodes);
   //          	cout << "thread sent nodes to master" << endl;
   //              payload.payloadStatus = Payload::WORKER_SHARED_NODES;
   //          	// share cuts to master. TODO: instead of move, keep a copy of local cuts.
   //          	// if (fst) payload.fCuts = move(fst.value());
   //          	// if (snd) payload.oCuts = move(snd.value());
   //          }
			// omp_unset_lock(&payload.lock);
		}
	}

	// display cut statistics. append to global string instead of printing.
	std::string msg = "Worker: " + std::to_string(id) + " processed " + std::to_string(nProcessed) + " nodes";
	cout << msg << endl;
}

void Inavap::DDSolver::Worker::shareCutsWithMaster(NodeExplorer &explorer, Inavap::DDSolver::Payload &payload) {
	// lock the paylaod and share cuts.
	// should share all the cuts to the payload.
	// payload.fCuts = move(explorer.feasibilityCuts);
	// if (explorer.optimalityCuts.size() > )
}


/**
 * Returns vector of nodes to process from the Payload.
 * @param status
 * @return
 */
vector<Inavap::Node> Inavap::DDSolver::Payload::getNodes(uint8_t &status) {

	/* Acquire lock and get already existing nodes (if any) from the payload. This is possible when
	 * the worker recently shared nodes with master, and master yet to take the nodes shared by this worker.
	 */

	if (omp_get_thread_num() ==1 && false) {

		cout << "Thread 1 reached here" << endl;
		if (omp_test_lock(&lock)) {
			cout << "test lock returned true" << endl;
			omp_unset_lock(&lock);
		}
		// there's a bug in the code if the test returns false, it should'nt happen.
		else cout << "test lock returned false" <<endl;
		// omp_unset_lock(&lock);
	}

	omp_set_lock(&this->lock);
	if (omp_get_thread_num() ==1)cout << "Thread 1 also reached here" << endl;
	if (!payloadNodes.empty()) {
		auto temp = std::move(payloadNodes);
		payloadStatus = WORKER_WORKING;
		omp_unset_lock(&this->lock);
		return temp;
	}
	payloadStatus = WORKER_NEEDS_NODES;

	int id = omp_get_thread_num();

	// cout << "Thread" << this_thread::get_id() << " Node Queue is empty. Indicating Master. " << endl;
	std::string message = "Thread " + std::to_string(id) + ": Queue is empty. Indicating Master."; cout << message << endl;
	// for now, set the status to finished
	// status = status | 1; return {};

	/* since openmp doesn't have conditional variables, I will create a spin lock manually.
	 * At this point, this thread still holds the lock.
	 */
	uint8_t st;

	int iteration = 0;

	// if (omp_get_thread_num() == 1) cout << "Iteration: ";
	while (true) {
		// sleep
		if (omp_get_thread_num() == 1) {
			// cout << iteration++ << " ";
		}

		st = this->payloadStatus;
		// omp_unset_lock(&this->lock); // if status changed?
		if (st == MASTER_ASSIGNED_NODES ||
			st == SOLVER_FINISHED)
			break;
		omp_unset_lock(&this->lock);;
		sleep(2);
		omp_set_lock(&this->lock);
	}
	// still holding lock
	if (st == SOLVER_FINISHED) {
		status = 1;
		omp_unset_lock(&this->lock);
		return {};
	}
	auto temp = std::move(payloadNodes);
	auto msg = "Thread " + std::to_string(id) + ": received " + std::to_string(temp.size()) + " nodes"; cout << msg << endl;
	status = 0;
	payloadStatus = WORKER_WORKING;
	omp_unset_lock(&this->lock);
	return temp;
}

void Inavap::DDSolver::startSolver() {

	cout << "Starting DD Solver." << endl;

	// create initial restricted tree and get cutset with desired max width.

	RestrictedDDNew restrictedDD{networkPtr, 128};
	Node root;
	auto cutset = restrictedDD.buildTree(root);

	if (cutset) cout << "Number of nodes from the first tree: " << cutset.value().size() << endl;

	vector<Node> cutsetNodes = cutset.value();
	// assume nodes is a vector of nodes.
	// divide the nodes to all workers.
	uint current = 0;
	// payloads = vector<Payload>(N_WORKERS);
	for (int x = 0; x < N_WORKERS; x++) {
		Payload p;
		p.id = x; // lock is initialized
		payloads.push_back(p);
	}
	for (const auto& node : cutsetNodes) {
		// round robin.
		auto index = (current++)%N_WORKERS;
		// omp_init_lock(&payloads[index].lock);
		payloads[index].payloadNodes.emplace_back(node);
	}

	// shared variables.
	#pragma omp parallel num_threads(N_WORKERS+1)
	{
		// #pragma omp single
		{
			// for (uint i = 0; i < N_WORKERS; i++) {
			//
			// 	// get paylaod init loc
			// 	omp_init_lock(&payloads[i].lock);
			// 	Worker worker {i, networkPtr};
			// 	workersGroup.push_back(worker);
			// 	// #pragma omp task untied
			// 	// worker.startWorker(this);
			// }
			//
			// cout << "all workers started" << endl;
			// // start master
			// Master m {networkPtr};
			// #pragma omp task untied
			// m.startMaster(*this);
			// cout << "Master finished" << endl;
		}

		#pragma omp single
		{
			// initialize payloads here.
			for (uint j = 0; j < N_WORKERS; j++) {
				omp_init_lock(&payloads[j].lock);
				Worker worker {j, networkPtr};
				workersGroup.push_back(worker);
			}
		}

		int thread_id = omp_get_thread_num();
		if (thread_id == N_WORKERS)
		{
			Master m {networkPtr};
			m.startMaster(*this);
			cout << "Master Finished" << endl;
		}
		else
			workersGroup[thread_id].startWorker(this);
	}
// 	for (uint i = 0; i < N_WORKERS; i++) {
// #pragma omp task untied
// 		{
// 			Worker worker {i, networkPtr};
// 			workersGroup.push_back(worker);
// 			worker.startWorker(this);
// 		}
// 	}
// 	// start master
// 	Master m{networkPtr};
// 	m.startMaster(std::ref(*this));
// 	// std::thread master (&Master::startMaster, &m, std::ref(*this));
// 	// if (master.joinable()) master.join();
// 	/* LATER: instead of creating new thread for master, let the main thread be the master thread */
// 	for (unsigned int i = 0; i < N_WORKERS; i++) {
// 		if (workers[i].joinable()) workers[i].join();
// 	}

	// destruct global cuts?
	cout << "Optimal Solution: " << optimal << endl;
	std::cout << std::flush;

}

void Inavap::DDSolver::NodeQueue::pushNode(Node node) {
	// push to local queue.
	pq.push(std::move(node));
}

void Inavap::DDSolver::NodeQueue::pushNodes(vector<Node> nodes) {
	for (auto&& node : nodes) pq.push(std::move(node));
}

Inavap::Node Inavap::DDSolver::NodeQueue::popNode() {
	// auto node = pq.top();
	Node node {pq.top()};
	pq.pop();
	return node;
}

vector<Inavap::Node> Inavap::DDSolver::NodeQueue::popNodes(size_t n) {
	vector<Inavap::Node> nodes;
	n = min(n, pq.size());
	nodes.reserve(n);
	while (n--) { nodes.emplace_back(pq.top()); pq.pop(); }
	return nodes;
}
