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
		auto [solution,execution_time] = solver.start(oracle);
		if (abs(solution-oracle) > 1e-5) cout << "ERROR: incorrect solution" << endl;
		cout <<"\n\n\n" << endl;
	}
}

void test_multiple(const std::string& filename, uint16_t threads, uint16_t n) {
	cout << "File: " << filename << endl;
	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{filename})};

	double oracle = solveStochasticModel(networkPtr.get());
	vector<double> exec_times(n);

	for (uint16_t i = 0; i < n; i++) {
		const auto now = std::chrono::system_clock::now();
		const auto t_c = std::chrono::system_clock::to_time_t(now);
		Inavap::DDSolver solver{networkPtr, threads};
		auto [solution, exec_time] = solver.start(oracle-10);
		if (abs(solution-oracle) > 1e-5) {
			cout << "ERROR: Incorrect solution" << endl.
			cout << "Aborting....";
			exit(-1);
		}
		exec_times[i] = exec_time;
		cout << "\n\n\n" << endl;
	}
	//
	double average_time = stats_mean(exec_times.data(), n);
	cout << "Mean execution time: " << average_time << endl;
}

int main(int argc, char* argv[]) {
	cout << "C++ version: " << __cplusplus << endl;

	string fileName;
	fileName = "/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/100_201_1_9.txt";
	// if (argc == 2) fileName = string(argv[1]);
	// fileName ="/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/40_93_20_10.txt";
	// test_single(fileName);
	// return 0;
	Network network{fileName};
	double optimal = solveStochasticModel(&network);
	// cout << "Optimal from stochastic model: " << optimal << endl;

	cout << "Reusing gurobi model between subproblems" << endl;
	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	Inavap::DDSolver solver{networkPtr, 8};
	auto [solution,exec_time ]= solver.start(optimal-10);
	return 0;
}

