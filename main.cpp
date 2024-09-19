#include <iostream>
#include "Network.h"
#include "grb.h"
#include "DD.h"

using namespace std;

void printNetwork(const Network& network){

	cout << "Number of nodes in the network: " << network.n << "\n"
		<< "Number of edges in the network: " << network.edges << "\n"
		<< "Number of nodes in v_bar: " << network.Vbar.size() << "\n";

	const auto& networkArc = network.networkArcs[2];

	cout << networkArc.arcId << " " << networkArc.tailId << " " << networkArc.headId << " "<< networkArc.rewards.size() << endl;

	/*for (const auto& e: network.Vbar)
		cout << e << " ";
	cout << endl; */

	//testWorking();

}

int main() {

	/* read input and build the core datastructures */
	/* assuming the input is text file with first line containing two numbers n (number of nodes) and number of edges */
	cout << "Starting program" << endl;
	string fileName = "/home/nandgate/Downloads/40_50_1.txt";
	cout << fileName << endl;
	Network network (fileName);

	// just for now.
	// remove 25 and add it at front
	network.Vbar.erase(std::remove(network.Vbar.begin(), network.Vbar.end(), 25), network.Vbar.end());
	network.Vbar.insert(network.Vbar.begin(), 25);

	//printNetwork(network);
	RestrictedDD dd{16};
	dd.build(network);

	cout << "Program completed." << endl;
}