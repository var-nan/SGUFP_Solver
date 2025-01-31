#include <iostream>
#include "Network.h"
#include "grb.h"
#include "DD.h"
#include <chrono>
#include "OriginalProblem.h"

using namespace std;

#include "DDSolver.h"

int main() {
	cout << "C++ version: " << __cplusplus << endl;
	string fileName ="C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_93_20_4.txt";
	Network network{fileName};

	SolveOriginalProblem(network);

	cout << "Solved Original problem " <<endl;
	cout << endl;
	cout << endl;

	const shared_ptr<Network> networkPtr{make_shared<Network>(Network{fileName})};
	cout << "Vbar order: "; for (auto id : networkPtr->Vbar) cout << id << " "; cout << endl;
	cout << "Max Width : " << MAX_WIDTH << endl;
	cout << "DEBUG enabled, solving for a subset of V Bar nodes. Single scenario." << endl;
	using std::chrono::high_resolution_clock;
	using std::chrono::duration_cast;
	using std::chrono::duration;
	using std::chrono::milliseconds;
	using std::chrono::seconds;
	auto t1 = high_resolution_clock::now();
	const auto now = std::chrono::system_clock::now();
	const auto t_c = std::chrono::system_clock::to_time_t(now);
	cout << endl << "Starting solver at " << std::ctime(&t_c);
	// DDSolver solver{networkPtr};
	Inavap::DDSolver solver{networkPtr, 1};
	// solver.startSolver();
	// solver.initialize();

//	int n_initial_cuts = 25;
//	auto cuts = solver.initializeCuts2(n_initial_cuts);
//	cout << "Number of initial cuts: " << n_initial_cuts << ". Optimality: " << cuts.second.cuts.size() <<
//		" , Feasibility: " << cuts.first.cuts.size() << endl;
	cout << "**********************************************************************************************************\n\n\n" << endl;

//	solver.startSolve(cuts);
	solver.startSolver();
	// solver.startPThreadSolver();

	auto t2 = high_resolution_clock::now();
	// cout << "Node queue strategy: LIFO" << endl;
	auto ms_int = duration_cast<seconds>(t2-t1);
	duration<double> ms_double = t2-t1;
	std::cout <<"program took " << ms_int.count() << " seconds" << endl;

	cout << "Solver finished" << endl;
	return 0;
}

