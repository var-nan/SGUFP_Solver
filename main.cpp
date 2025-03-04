#include <iostream>
#include "Network.h"
#include "grb.h"
#include "DD.h"
#include <chrono>
#include "OriginalProblem.h"

using namespace std;

#include "DDSolver.h"

void test_single(const std::string& fileName) {
	cout << "Testing " << fileName << endl;
	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};

	SolveOriginalProblem(*networkPtr);

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
		solver.startSolver();
		auto t2 = high_resolution_clock::now();
		auto ms_int = duration_cast<seconds>(t2-t1);
		duration<double> ms_double = t2-t1;
		std::cout <<"Solver finished: " << ms_int.count() << " seconds" << endl;
		solver.printWorkerStats();
	}
}

int main() {

	// std::vector<std::string> files = {
	//
	// };
	//
	// return 0;
	cout << "C++ version: " << __cplusplus << endl;
	string fileName ="/mnt/c/Users/nandgate/CLionProjects/SGUFP_Solver/40_91_20_1.txt";
	Network network{fileName};

	SolveOriginalProblem(network);

	cout << "Solved Original problem " <<endl;
	cout << endl;
	cout << endl;

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
	solver.startSolver();
	auto t2 = high_resolution_clock::now();
	auto ms_int = duration_cast<seconds>(t2-t1);
	duration<double> ms_double = t2-t1;
	std::cout <<"Solver finished: " << ms_int.count() << " seconds" << endl;
	solver.printWorkerStats();
	return 0;
}

