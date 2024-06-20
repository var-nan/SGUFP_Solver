//
// Created by nandgate on 6/3/24.
//

#include "Solver.h"

/**
 * Creates W_bar data structure from the given input path.
 * @param path: vector<uint32_t>, contains arc decisions that maximizes the objective from root to terminal.
 * @param network: current network.
 *
 * @return: W_bar:
 */
static inline vector<vector<uint32_t>> convertPathToW(const Network& network, vector<uint32_t>& path){
	// converts given path from DD to W data structure. path is a vector of arc labels from root to terminal.
	// order of path. root to terminal. size of path vector is size of tree-1
	const auto& networkNodes = network.networkNodes;
	//const auto& networkArcs = network.networkArcs;
	const auto& Vbar = network.Vbar;

	vector<vector<uint32_t>> w_bar; // TODO: initialize w_bar.

	uint32_t last_accessed = 0;
	for (uint32_t q = 0; q < network.Vbar.size(); q++){
		const auto q_node = networkNodes[Vbar[q]];
		for (uint32_t i = 0; i < q_node.inDegree; i++){
			//const auto arc = q_node.inNodeIds[i];
			const auto incomingNode = q_node.inNodeIds[i];
			w_bar[q][i] = path[last_accessed++];
		}
	}
	return w_bar;
}

/**
 * converts y variable to w variable.
 */
inline void convertToW(const vector<vector<vector<bool>>>& y, vector<vector<uint32_t>>& w){
	// w is a 2 dimensional vector of values (q,h), y is a 3D vector of values (q, i,j)
	// w[Vbar.size(), indexset.size()]. first two dimensions should be same for w and y.
	// TODO: INITIALIZE y and w variables with valid dimensions.
	for (uint32_t q = 0; q < w.size(); q++){
		for (uint32_t h = 0; h < w[q].size(); h++){
			int t = -1;
			for (uint32_t j = 0; j < y[q][h].size(); j++){
				if (y[q][h][j]) {
					t = j;
					break;
				}
			}
			if (t == -1) t = 0;
			w[q][h] = t;
		}
	}
}

inline int sign(int x){
	if (x > 0) return 1;
	else if (x == 0) return 0;
	else return -1;
}


inline void convertToY(const Network& network, vector<vector<uint32_t>>& w, const vector<vector<vector<bool>>>& y) {
	// y^q_ij = 1 - sign (|w^q_ij - ind(q,j)|)

	for (uint32_t q_index = 0; q_index < y.size(); q_index++) {
		const auto q = network.Vbar[q_index];
		const auto q_node = network.networkNodes[q];

		for (uint32_t i = 0; i < y[q].size(); i++) {
			for (uint32_t j = 0; j < y[q][i].size(); j++) {

			}

		}
	}
}

void Solver::solveStochasticDD() {
	/* start of algorithm described in the paper */

	// initialize set of partial assignments

	/*
	 * Build a restricted DD with global order among the variables.
	 * get longest path from terminal to root in w-variable space.
	 * transform w-variables to y-space.
	 * solve Bender's sub-problem using y to obtain benders cut C.
	 * If C is new cut, add it to list of cuts and refine the DD with C.
	 *
	 */


	uint32_t  w_opt = 0;
	uint32_t w_current = 0;

	const auto& network = this->network;

	GRBEnv environment;
	GurobiSolver solver{environment, network.n};

	DDNode root;
	uint32_t maxWidth = 64;
	queue<DDNode*> nodeQueue;

	// scope for parallelization.
	while (!nodeQueue.empty()) {
		auto current = *(nodeQueue.front());
		nodeQueue.pop();

		// build restricted DD.
		RestrictedDD restrictedDD {maxWidth};
		restrictedDD.build(network, current.nodeId);
		int t = 0;
		while (t < 4) {
			t++;
			// get longest path from restricted DD.
			vector<uint32_t> path;
			path = restrictedDD.getPath();
			vector<vector<uint32_t>> w = convertPathToW(network, path);

			// TODO: code to compute longest path and length.

			// convert w to Y.
			vector<vector<vector<bool>>> y_bar;

			// solve subproblem using benders.
			for (uint32_t scenario = 0; scenario < network.nScenarios; scenario++){
				solver.solveSubProblem(network, y_bar, scenario); // either returns optimality cut or feasibility cut.
			}

			// ASAP define Cut datastructure and get cut.



			// refine DD wrt. to cut.

		}

		// if w_bar is greater than w_opt, then increase the w_opt
		if (w_current > w_opt){
			w_opt = w_current;
			// update y_opt and z_opt
		}

		// if restricted DD is exact, then skip
		if (!restrictedDD.isExact){
			// build relaxed DD
			RelaxedDD relaxedDD;
			relaxedDD.build(network, current.nodeId);
			// find longest path and
			uint32_t w_current_relaxed = 10;

			if (w_current_relaxed > w_opt){
				// solve subproblem using y_bar to obtain bender's cut.
				//

				// add exact cutset to queue.
				for (auto& node: relaxedDD.getCutSet()){
					nodeQueue.push(&node);
				}
			}
		}
	}
}
