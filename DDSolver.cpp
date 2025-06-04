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

	if (nfCuts > FEASIBILITY_CONTAINER_CAPACITY) {
		uint nContainers = static_cast<uint> (ceil((1.0*nfCuts)/FEASIBILITY_CONTAINER_CAPACITY));
		for (auto i = 0; i < nContainers; i++)
			feasibilityPointers.push_back(new CutContainer(FEASIBILITY_CONTAINER_CAPACITY));

        uint c_index = 0;
        uint n_so_far = 0;
        for (auto& container : cuts.first) {
            for (auto& cut : container) {
                feasibilityPointers[c_index]->insertCut(cut);
            	if (!(++n_so_far % FEASIBILITY_CONTAINER_CAPACITY)) c_index++;
            }
        }
	}

	size_t nOCuts = accumulate(cuts.second.begin(), cuts.second.end(), 0,
		[](size_t sum, const CutContainer &cont) { return sum + cont.size(); });
	if (nOCuts > OPTIMALITY_CONTAINER_CAPACITY) {
		uint nContainers = static_cast<uint> (ceil( static_cast<double>(nOCuts)/OPTIMALITY_CONTAINER_CAPACITY));
		for (auto i = 0; i < nContainers; i++)
			optimalityPointers.push_back(new CutContainer(OPTIMALITY_CONTAINER_CAPACITY));
        uint c_index = 0;
        uint n_so_far = 0;
        for (auto& container : cuts.second) {
            for (auto& cut: container) {
                optimalityPointers[c_index]->insertCut(cut);
            	if (!(++n_so_far % OPTIMALITY_CONTAINER_CAPACITY)) c_index++;
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
	// // worker-mode.
	// double zOpt = solver.optimal.load(memory_order::relaxed);
	// size_t nProcessed = 0;
	//
	// while (!nodeQueue.empty() && n--) {
	// 	auto result = explorer.process(nodeQueue.popNode(), zOpt, solver.feasCutsGlobal, solver.optCutsGlobal);
	// 	nProcessed++;
	// 	// TODO: what to do with local cuts?
	// 	if (result.status == OutObject::STATUS_OP::SUCCESS) {
	// 		if (result.lb > zOpt) {
	// 			while ((!solver.optimal.compare_exchange_weak(zOpt, result.lb, memory_order::relaxed)) && (zOpt < result.lb));
	// 		}
	// 		if (result.ub > zOpt && !result.nodes.empty()) nodeQueue.pushNodes(result.nodes);
	// 	}
	// }
	// return nProcessed;
	return 0;
}

void Inavap::DDSolver::Master::startMaster(DDSolver &solver) {

	NodeExplorer explorer{networkPtr};
	size_t nProcessed = 0;

	const uint16_t N_WORKERS = solver.N_WORKERS;

	vector<Payload::STATUS> worker_status(N_WORKERS, Payload::WORKER_WORKING); // 0 - sleeping, 1- working
	// maintain payload status.
	for (;;) {
		// number of idle workers, and number of processing workers.
		uint16_t idle = 0, processing = 0;

		// iterate through all workers.
		for (uint16_t i = 0; i < N_WORKERS; i++) {
			auto& worker = solver.payloads[i];
			auto st = worker.payloadStatus.load(std::memory_order::acquire);

			if (st == Payload::WORKER_NEEDS_NODES) { // worker should be in waiting state.
				// check local queue and send them to queue if queue is not empty
				// cout<< "master: Worker " << i << " needs nodes" << endl;
				size_t q_size = nodeQueue.get_size();
				if (q_size) {
					size_t k = (q_size &1) ? (1+ q_size>>1) : q_size>>1; // TODO; fill this.
					cout << "sharing " << k << " nodes to worker" << endl;
					// assign nodes to worker.
					llist new_list = nodeQueue.pop_k(k);

					{
						size_t x = 0;
						for (const lf_node *begin = new_list.start; (begin) ; begin = begin->next) x++;
						assert(x == new_list.n && x > 0);
						if (!new_list.start || !new_list.end || !new_list.n) {cout << "Empty list" << endl;}
					}

					worker.private_queue.push(new_list);
					cout << "Master shared nodes" << endl;

					// signal worker
					worker.payloadStatus.store(Payload::WORKER_WORKING);
					worker.payloadStatus.notify_one();
					worker_status[i] = Payload::WORKER_WORKING;
					cout << "master signaled worker" << endl;
					processing++;
				}
				else {
					// cout << "master queue is empty" << endl;
					idle++;
					worker_status[i] = Payload::WORKER_NEEDS_NODES;
				}
			}
			else {
				// status is either WORKER_WORKING or MASTER_ASSIGNED_NODES
				worker_status[i] = Payload::WORKER_WORKING;
				processing++;
			}
		}

		if (idle == N_WORKERS && nodeQueue.empty()) {
			// solver is finished, all the threads must be waiting now.
			cout << "master: solver is finished, indicating workers." << endl;
			for (uint16_t i = 0; i < N_WORKERS; i++) {
				auto& worker = solver.payloads[i];
				worker.payloadStatus.store(Payload::SOLVER_FINISHED);
				// std::atomic_notify_one(&worker.payloadStatus);
				worker.payloadStatus.notify_one();
			}
			break;
		}

		if (!idle && !nodeQueue.empty()) {
			// process nodes
		}


		for (uint16_t i = 0; idle && i < N_WORKERS; i++) {
			if (worker_status[i] == Payload::WORKER_NEEDS_NODES) continue;
			// take nodes from this worker.
			auto& worker = solver.payloads[i];
			auto st = worker.payloadStatus.load(memory_order::acquire); // this is not needed.
			if (st == Payload::WORKER_NEEDS_NODES) continue; // this is not needed
			llist new_list = worker.private_queue.m_pop(0.4);
			if (!new_list.start) continue; // either this worker is empty or not have enough nodes.
			// push to local queue.
			cout << "fetched " << new_list.n << " nodes from " << i << endl;
			// check if new_list  has n nodes or not.
			{
				size_t x = 0;
				for (const lf_node *begin = new_list.start; (begin); begin = begin->next) x++;
				if (x != new_list.n) {
					cout << "values mismatch: " << endl;
					exit(-1);
				}
			}
			nodeQueue.push(new_list);
		}
	}
	// post-completion tasks by master?
	// solver.CutResources.printStatistics();
}

void Inavap::DDSolver::Worker::startWorker(DDSolver *solver) {
	cout << "Worker started" << endl;

	/* function: returns true if it is time to read the atomic variable of global optimal, else increases counter. */
	auto is_poll_time = [](size_t &counter) { // resets counter when reaches limit.
       if (counter == POLL_FREQUENCY) { counter = 0; return true; }
        counter++; return false;
    };

	auto is_check_time = [](uint &counter) {
		if (counter == PAYLOAD_CHECK_TICKS) {
			counter = 0;
			return true;
		}
		counter++;
		return false;
	};

	NodeExplorer explorer{networkPtr};
	auto& payload = solver->payloads[id];
	payload.id = id;
	auto& private_queue = payload.private_queue;
	double zOpt = solver->optimal.load(memory_order::relaxed);

	uint8_t done = 0; // flag to indicate solver is finished.
	// size_t nProcessed = 0; // # nodes processed so far.
	size_t counter = 0;
	uint nOptShared = 0;		// # local optimality cuts shared with master so far.
	uint nFeasShared = 0;		// # local feasibility cuts shared with master so far.
	uint prevLocalCount = 0;	//
	uint prevCount = 0;			// # newly generated cuts since last shared with global.
	uint payloadCheckTick = 0;	//

	/* Operations on same variable by a single thread with relaxed-order will obey the happens-before relationships.
	 * The access to a single atomic variable from  the same thread can't be reordered. Once a given thread see a
	 * particular value of atomic variable, any subsequent read by that thread is guaranteed to retrieve same value or
	 * latest value since the last read operation.
	 */

	/* Still don't know how the atomic store is implemented at the hardware level. Does it directly write to memory ?
	 * or does it store in caches in case of relaxed ordering? I think I know bit better now.
	 */
	while (private_queue.empty()){} // after this point, the worker should see a valid head pointer.
	cout << "Worker seen head pointer" << endl;

	for (;;) {

		for (lf_node *node_ptr; (node_ptr = private_queue.pop()); ) {

			Node node = *node_ptr;
			// auto str = "node: " + to_string(node.globalLayer); cout << str << endl;
			if (node.ub <= zOpt) {
				nPrunedByBound++;
				delete node_ptr; // should not cause leak.
				continue;
			}
			auto result = explorer.process(node, zOpt, solver->feasCutsGlobal, solver->optCutsGlobal);

			delete node_ptr; // free nodes.
			#ifdef SOLVER_COUNTERS
			nProcessed			+= (result.status == OutObject::STATUS_OP::SUCCESS);
			nFeasibilityPruned	+= (result.status == OutObject::STATUS_OP::PRUNED_BY_FEASIBILITY_CUT);
			nOptimalityPruned	+= (result.status == OutObject::STATUS_OP::PRUNED_BY_OPTIMALITY_CUT);
			#endif
			nQueue				+= (result.nodes.size());
			if (result.status == OutObject::STATUS_OP::SUCCESS) {
				//
				if (result.lb > zOpt) {
					/* short-circuit operation. If the store is successful, exit the loop, else keep trying until
					 * either one of holds (1) store successful, (2) another best optimal found.
					 * Can be replaced with compare_exchange_strong, but might need additional check. */
					while (!solver->optimal.compare_exchange_weak(zOpt, result.lb, memory_order::relaxed)
						&& (zOpt < result.lb)) {}

					// zOpt is not updated after CAS operation.
					if (solver->optimal.load(memory_order::acquire) == result.lb) zOpt = result.lb;
#ifdef SOLVER_COUNTERS
					 if (result.lb == zOpt) {
					 	const auto now = std::chrono::system_clock::now();
					 	const auto t_c = std::chrono::system_clock::to_time_t(now);
					 	cout << "Thread: " << id << " , optimal lb: " << zOpt << " set at, " << std::ctime(&t_c) << endl;
					 }
#endif
					counter = 0;
				}
				else if (is_poll_time(counter)) { // check periodically for latest optimal value.
					zOpt = solver->optimal.load(memory_order::relaxed);
				}
				if (result.ub > zOpt && !result.nodes.empty()) {
					llist nodes = convert(result.nodes);
					assert((nodes.start && nodes.end && nodes.n));
					private_queue.push(nodes);
				}
			}
		}

		#ifdef SOLVER_COUNTERS
			sleepTimes++;
			const auto start = std::chrono::system_clock::now();
		#endif

		cout << "worker needs nodes"<< endl;
		// private queue is empty. indicate to master.
		payload.payloadStatus.store(Payload::WORKER_NEEDS_NODES, std::memory_order::release);

		// sleep until master signals.
		payload.payloadStatus.wait(Payload::WORKER_NEEDS_NODES);
		cout <<"Woke up by master" << endl;

		// std::atomic_thread_fence(std::memory_order::acquire);
		#ifdef SOLVER_COUNTERS
			const auto end = std::chrono::system_clock::now();
			sleepDuration += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
		#endif

		// either solver is finished or master assigned some nodes.
		auto status = payload.payloadStatus.load(std::memory_order::acquire);
		// while (status == Payload::WORKER_NEEDS_NODES) {status = payload.payloadStatus.load(memory_order::acquire);}
		if (status == Payload::SOLVER_FINISHED) {
			break;
		}
	}

	cout << "finished: " << payload.id << endl;
}

/**
 * Starts the solver with the given value as the current optimal.
 * @param known_optimal
 */
double Inavap::DDSolver::startSolver(double known_optimal) {

	optimal.store(known_optimal, memory_order::relaxed);

	// create initial restricted tree and get cutset with desired max width.

	RelaxedDDNew relaxedDD{networkPtr.get()};
	Node root;
	relaxedDD.buildTree(root);
	auto cutset = relaxedDD.getCutset(DOUBLE_MAX);

	vector<Node> cutsetNodes = cutset;
	uint current = 0;
	payloads = vector<Payload>(N_WORKERS);

	std::reverse(cutsetNodes.begin(), cutsetNodes.end());
	llist nodes = convert(cutsetNodes);
	lf_node *st = nodes.start, *end = nodes.end;
	size_t sz = nodes.n;
	for (auto& payload : payloads) {
		llist worker_nodes = {st, st, 1};
		lf_node *temp = st->next;
		payload.private_queue.push(worker_nodes);
		st = temp;
		assert(st != nullptr);
		sz--;
	}
	// check if nodes have correct nodes
	size_t finalsize = 0;
	for (const lf_node *begin = st; (begin); begin = begin->next){finalsize++;}
	if (finalsize != sz) {
		cout << "values mismatch"<< endl;
		abort();
	}
	nodes = {st, end, sz};

	// payloads are ready. launch worker threads.
	for (uint i = 0; i < N_WORKERS; i++) {
		Worker worker{i, networkPtr};
		workersGroup.push_back(worker);
		workers.emplace_back(&Worker::startWorker, &workersGroup[i], this);
	}
	// let main thread work as master, instead of spawning new thread.
	Master m{networkPtr, nodes};
	m.startMaster(std::ref(*this));

	for (unsigned int i = 0; i < N_WORKERS; i++) {
		if (workers[i].joinable()) workers[i].join();
	}

	return optimal.load();

	//  TODO: post processing needed. (clean up resources).

}

double Inavap::DDSolver::start(double known_opt) {

	const auto startTime = std::chrono::high_resolution_clock::now();
	double solution = startSolver(known_opt);
	const auto endTime = std::chrono::high_resolution_clock::now();
	auto duration_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);

#ifdef SOLVER_COUNTERS
	printWorkerStats();
#endif

	// compute total nodes.
	size_t totalQueueNodes = 0;
	for (const auto& worker : workersGroup) totalQueueNodes += worker.nQueue;

	std::cout << "Optimal solution: " << solution
		<< ". Explored "<< totalQueueNodes << " nodes (entire search space) in " << duration_seconds.count() << " seconds." << endl;
	std::cout << "Threads : " << N_WORKERS << " Workers and 1 Master." << std::endl;
	return solution;
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
