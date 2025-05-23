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

	double oracle = solveStochasticModel(networkPtr.get());

	for (uint16_t i = 1; i <= 16; i*= 2) {
		const auto now = std::chrono::system_clock::now();
		const auto t_c = std::chrono::system_clock::to_time_t(now);
		cout << "Starting solver " << std::ctime(&t_c) << endl;
		Inavap::DDSolver solver{networkPtr, i};
		double solution = solver.start(oracle);
		if (abs(solution-oracle) > 1e-5) cout << "ERROR: incorrect solution" << endl;
		cout <<"\n\n\n" << endl;
	}
}

int main(int argc, char* argv[]) {
	cout << "C++ version: " << __cplusplus << endl;

	string fileName;
	fileName = "/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/40_87_20_3.txt";
	// if (argc == 2) fileName = string(argv[1]);
	// fileName ="/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/40_93_20_10.txt";
	// test_single(fileName);
	// return 0;
	Network network{fileName};
	double optimal = solveStochasticModel(&network);
	cout << "Optimal from stochastic model: " << optimal << endl;

	cout << "Reusing gurobi model between subproblems" << endl;
	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	Inavap::DDSolver solver{networkPtr, 4};
	double solution = solver.start(optimal-100000);
	return 0;
}

