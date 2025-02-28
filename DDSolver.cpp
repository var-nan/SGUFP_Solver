//
// Created by nandgate on 10/24/2024.
//

#include "DDSolver.h"
#include <random>
#include <chrono>

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

void DDSolver::initialize(double opt) {
    // place root node to queue
    Node_t node;

    node.lb = std::numeric_limits<int64_t>::min();
	node.ub = 5000000;
    node.globalLayer = 0;
	setLB(opt);
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
	int counter = 0;
    while (!nodeQueue.empty()) { // conditional wait in parallel version

        Node_t node = nodeQueue.getNode();
        // cout << "Procesisng Node from layer: "<< node.globalLayer << " LB: " << node.lb << " , UB: " << node.ub << " global: " << getOptimalLB()<< endl;
        #ifdef DEBUG
        // cout << "Processing node from layer: " << node.globalLayer << " lb: " << node.lb << " , ub: " << node.ub;
        cout << " . global lower bound: " << getOptimalLB() << endl;
        #endif
        if (node.ub <= getOptimalLB()) {
            #ifdef SOLVER_STATS
            numPrunedByBound++;
            // cout << "Pruned by bound." << endl;
            #endif
            continue; // look for another
        }

        // start node processor
    	// cout << "process 4 stared"<< endl;
    	OutObject result = {0,0,{},0};

    	// if (counter >= node.globalLayer ) {
    	// 	counter =0;
    	// 	// cout << "_4_";
    	// 	result = explorer.process4(node, (getOptimalLB() )); // use co-routines to update globalLB in between.
    	// }else {
    	// 	counter += 1;
    	// 	// cout << "_5_";
    	// 	result = explorer.process5(node, (getOptimalLB() )); // use co-routines to update globalLB in between.
    	// }
    	// cout << "zart" << endl;
    	result = explorer.processX4(node, getOptimalLB() ); // use co-routines to update globalLB in between.
    	cout << "new UB:  " << result.ub << "   ";
    	cout << result.nodes.size() << endl;
    	// cout << "process 4 finished"<< endl;
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
    return optimalLB;
}

void DDSolver::setLB(double lb) {
    optimalLB = lb;
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
	// cout << "zart1" << endl;
	// auto vecc = networkPtr->Vbar;
	// for(auto item:vecc) {
	// 	cout << item << "*" ;
	// }
	// cout << endl;
	// std::reverse(vecc.begin(), vecc.end());
	// for(auto item:vecc) {
	// 	cout << item << "*" ;
	// }
	// cout << endl;
	// networkPtr->Vbar = vecc;
	// cout << "zart1" << endl;
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
        // cout << "Solution selected: "; for (auto s: solution) cout << s <<" "; cout <<endl;
        solutions.push_back(solution);
        // cout << "Solution: "; for (auto sol: solution) cout << sol << " "; cout << endl;
        // get cut
        GuroSolver solver{networkPtr, env};
        const auto& y_bar = w2y(solution, networkPtr);
        auto cut = solver.solveSubProblemInstance(y_bar); // LATER set random scenario.
        if (cut.cutType == FEASIBILITY) {
            fCuts.insertCut(cut);
        	cout << "feas cut: " << endl;
        }
        else oCuts.insertCut(cut);
    		cout << "opt cut" << endl;
        n--;
    }
    return make_pair(fCuts, oCuts);
}

pair<CutContainer, CutContainer> DDSolver::initializeCuts3() {

    // start with root and find random paths to terminal.
    vector<vi> solutions;
    CutContainer fCuts{FEASIBILITY};
    CutContainer oCuts{OPTIMALITY};
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

    // solutions.reserve(n);
	for (auto item : networkPtr->Vbar) {
		cout << item << " * " ;
	}
	cout << endl;
	int ct = 0;
	reverse(networkPtr->Vbar.begin(), networkPtr->Vbar.end());
	for (auto item : networkPtr->Vbar) {
		cout << item << " * ";
	}
	cout << endl;
    while (ct < networkPtr->Vbar.size()) {
    	cout << "first node:" <<networkPtr->Vbar[0]<< endl;

    	ct += 1;
    	DD restrictedDD{networkPtr, RESTRICTED};
    	DDNode root{0};
    	restrictedDD.build(root);


    	auto lowerBound = std::numeric_limits<int64_t>::min();


    	double previousLowerBound;
    	vi previousSolution;
    	while (true) {
    		vi solution = restrictedDD.solution();
    		if (previousSolution == solution) {
    			if (previousLowerBound == lowerBound) {
    				break;
    			}
    			previousLowerBound=lowerBound;
    		}else {
    			previousSolution = solution;
    			previousLowerBound = lowerBound;
    		}
    		auto y_bar = w2y(solution , networkPtr);
    		GuroSolver solver{networkPtr, env};
    		auto cut = solver.solveSubProblemInstance(y_bar);
    		if (cut.cutType == FEASIBILITY) {
    			fCuts.insertCut(cut);
    			if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
    				// cout << "exact tree pruned by feas actual ref!" << endl;
    				break;
    			}
    		}else {
    			oCuts.insertCut(cut);
    			lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
    			if (lowerBound <= optimalLB) {
    				// cout << "exact tree pruned by opt actual ref!" << endl;
    				break;
    			}
    		}
    	}
    	networkPtr->vbarReOrder(1);
    	cout << "fCuts: " << fCuts.cuts.size() << endl;
    	cout << "oCuts: " << oCuts.cuts.size() << endl;
    }
	reverse(networkPtr->Vbar.begin(), networkPtr->Vbar.end());
	networkPtr->vbarReOrder(0);
	for (auto item : networkPtr->Vbar) {
		cout << item << " * ";
	}
	cout << endl;
    return make_pair(fCuts, oCuts);
}

pair<CutContainer, CutContainer> DDSolver::initializeCuts4() {

    // start with root and find random paths to terminal.
    vector<vi> solutions;
    CutContainer fCuts{FEASIBILITY};
    CutContainer oCuts{OPTIMALITY};
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_OutputFlag, 0);

	int ct = 0;
    while (ct == 0) {
    	ct += 1;
    	DD restrictedDD{networkPtr, RESTRICTED};
    	DDNode root{0};
    	restrictedDD.build(root);
    	auto lowerBound = std::numeric_limits<int64_t>::min();
    	double previousLowerBound;
    	vi previousSolution;
    	while (true) {
    		vi solution = restrictedDD.solution();
    		if (previousSolution == solution) {
    			if (previousLowerBound == lowerBound) {
    				break;
    			}
    			previousLowerBound=lowerBound;
    		}else {
    			previousSolution = solution;
    			previousLowerBound = lowerBound;
    		}
    		auto y_bar = w2y(solution , networkPtr);
    		GuroSolver solver{networkPtr, env};
    		auto cut = solver.solveSubProblemInstance(y_bar);
    		if (cut.cutType == FEASIBILITY) {
    			fCuts.insertCut(cut);
    			if (!restrictedDD.applyFeasibilityCutRestrictedLatest(cut)) {
    				break;
    			}
    		}else {
    			oCuts.insertCut(cut);
    			lowerBound = restrictedDD.applyOptimalityCutRestrictedLatest(cut);
    		}
    	}

    }
	// reverse(networkPtr->Vbar.begin(), networkPtr->Vbar.end());
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
            auto cut = solver.solveSubProblemInstance(y_bar);
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
    // vector<Node_t> t1, t2;

	{

		// single iteration of node explorer.
		Node_t node {{},{}, numeric_limits<double>::lowest(), numeric_limits<double>::max(),0};
		NodeExplorer explorer {networkPtr};
		cout << "process 4 start"<< endl;
		auto result = explorer.process4(node, numeric_limits<double>::lowest());
		cout << "process 4 finished"<< endl;
		if (result.lb > globalLB.load(memory_order_acquire))
			globalLB.store(result.lb, memory_order_release);
//		nodeQueue.pushNodes(result.nodes);
		int size = result.nodes.size();
		for (int i = 0; i < size; i++) {
		    if (i < size/2) {tempQ1.pushNode(result.nodes[i]); }//t1.push_back(result.nodes[i]);}
		    else {tempQ2.pushNode(result.nodes[i]);} //t2.push_back(result.nodes[i]);}
		}
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

    // vector<WorkerElement> workers(2);
    // workers[0].addNodes(t1); workers[1].addNodes(t2);
	std::thread worker{&DDSolver::processWork, this, tempQ1, cuts};
	std::thread worker2{&DDSolver::processWork, this, tempQ2, cuts};

	if (worker.joinable()) {
		//cout << "Worker thread is joinable" << endl;
		worker.join();
	}
	if (worker2.joinable()) worker2.join();
	//else {cout << "worker thread is not joinable" << endl;}

	//assert(nodeQueue.empty());
	//cout << "Number of nodes in the queue: " << nodeQueue.size() << endl;

	cout << "Optimal Solution : " << globalLB.load() << endl;

}


void DDSolver::processWork(NodeQueue q, pair<CutContainer, CutContainer> cuts) {
	cout << "thread: "<< this_thread::get_id << " starting" << endl;

	NodeExplorer explorer{networkPtr, cuts}; // get initial cuts later.
    //
	// int id = 0;
	// auto& payload = workers[id];
	// auto nodeVec  = payload.getNodes();
	NodeQueue localQueue = q; // initialize localQueue later.

	double zOpt = globalLB.load(memory_order_acquire);
	size_t nProcessed  = 0;

	// while (!isCompleted.load(memory_order_acquire)) {
		// if (localQueue.empty()) {
		// 	auto nodes = payload.getNodes();
		// 	localQueue.pushNodes(nodes);
		// }
		while (!localQueue.empty()){
			Node_t node = localQueue.getNode();
			// cout << "thread: " << this_thread::get_id <<"processing node: " << node.globalLayer << endl;
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

			// //if master wants some work?
			// if (payload.masterRequireNodes()) {
			// 	// share half of nodes with master.
			//
			// }

		}
	// }

	// display stats.
	cout << "thread: " << this_thread::get_id() << " processed " << nProcessed << " nodes." << endl;
	explorer.displayCutStats();
}

void Payload::addNodes(vector<Node_t> nodes) {
	{
		scoped_lock l{lock};
		nodes_ = move(nodes);
	}
	status.store(WORKER_WORKING, memory_order_release);
	cv.notify_one();
}

vector<Node_t> Payload::getNodes() {
	// return nodes from the master nodes
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
	status.store(WORKER_NEEDS_NODES,memory_order_release);
	while (true) {
		std::unique_lock ul{lock};
		cv.wait(ul, [&]{return !nodes_.empty();});
		vector<Node_t> work = std::move(nodes_);
		status.store(WORKER_WORKING, memory_order_release);
		// nodes_.clear();
		return work;
	}
}

bool Payload::masterRequireNodes() const noexcept {
	return status.load(memory_order_relaxed) & MASTER_NEEDS_NODES;
}

