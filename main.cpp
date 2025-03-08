#include <iostream>
#include "Network.h"
#include "grb.h"
#include "DD.h"
#include <chrono>
#include "OriginalProblem.h"

#include "include/StochasticModel.h"

using namespace std;

#include "DDSolver.h"

void test_single(const std::string& fileName) {
	cout << "Testing " << fileName << endl;
	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};

	double solution = solveStochasticModel(networkPtr.get());

	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	using std::chrono::seconds;

	for (uint16_t i = 1; i <= 16; i*= 2) {
		const auto now = std::chrono::system_clock::now();
		const auto t_c = std::chrono::system_clock::to_time_t(now);
		cout << endl << "Starting solver with " << i <<" threads at" << std::ctime(&t_c);
		auto t1 = high_resolution_clock::now();
		Inavap::DDSolver solver{networkPtr, i};
		solver.startSolver(solution);
		auto t2 = high_resolution_clock::now();
		auto ms_int = duration_cast<seconds>(t2-t1);
		duration<double> ms_double = t2-t1;
		std::cout <<"Solver finished: " << ms_int.count() << " seconds" << endl;
		solver.printWorkerStats();
	}
}

int main() {


	// vector<string> files = {};
	// for (const auto file : files) {
	// 	test_single(file);
	// }
	// return 0;
	cout << "C++ version: " << __cplusplus << endl;
	string fileName ="/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/40_93_20_2.txt";
	Network network{fileName};

	double optimal = solveStochasticModel(&network);
	cout << "Optimal from stochastic model: " << optimal << endl;


	/*
	{
		// compile new relaxed DD
		Inavap::RelaxedDDNew newDD{&network};
		Inavap::Node node{{},{},
			std::numeric_limits<double>::lowest(),
			std::numeric_limits<double>::max(),0};

		newDD.buildTree(node);

		cout << "Relaxed DD Nodes:  " << newDD.nodes.size() << endl;
		cout << "Relaxed DD arcs: " << newDD.arcs.size() << endl;
		cout << "Relaxed DD tree: " << newDD.tree.size() << endl;
		for (const auto& layer : newDD.tree) cout << layer.size() << " ";
		cout << endl;

		// print solution
		cout << "Solution: " << endl;
		auto sol = newDD.getSolution();
		for (auto s : sol) cout << s << " "; cout << endl;

		// compare it with the old relaxed DD
		Inavap::RelaxedDD relaxedDD{make_shared<Network>(network)};
		relaxedDD.buildTree(node);
		cout << "Old RelaxedDD Nodes: " << relaxedDD.nodes.size() << endl;
		cout << "arcs: " << relaxedDD.arcs.size() << endl;
		cout << "Tree: " << relaxedDD.tree.size() << endl;

		cout << "new Tree built" <<endl;
		return 0;
	}
	cout << "Solved Original problem " <<endl;
	cout << endl;
	cout << endl;

	*/

	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	using std::chrono::seconds;
	auto t1 = high_resolution_clock::now();
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	Inavap::DDSolver solver{networkPtr, 3};
	solver.startSolver(optimal);
	auto t2 = high_resolution_clock::now();
	auto ms_int = duration_cast<seconds>(t2-t1);
	duration<double> ms_double = t2-t1;
	std::cout <<"Solver finished: " << ms_int.count() << " seconds" << endl;
	solver.printWorkerStats();
	return 0;
}

