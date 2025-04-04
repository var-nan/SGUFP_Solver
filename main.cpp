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

int main(int argc, char* argv[]) {
	cout << "C++ version: " << __cplusplus << endl;

	string fileName;
	// fileName ="/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/instances/40_91_20_1.txt";
	if (argc == 2) fileName = string(argv[1]);
	Network network{fileName};
	double optimal = solveStochasticModel(&network);
	cout << "Optimal from stochastic model: " << optimal << endl;

	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	Inavap::DDSolver solver{networkPtr, 8};
	double solution = solver.start(optimal-10);
	return 0;
}

