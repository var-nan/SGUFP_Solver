//
// Created by nandgate on 1/6/2025.
//
// test cases for Inavap namespace.

#include "test_utils.h"

using namespace Inavap;

// class NodeExplorerTest: public testing::Test {
// 	// test for node explroer class.
// protected:
// 	// Inavap::NodeExplorer explorer;
// 	const string fileName {"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_91_20_1.txt"};
// 	const shared_ptr<Network> networkPtr = make_shared<Network>(Network{fileName});
// 	Inavap::NodeExplorer explorer{networkPtr};
// public:
// 	NodeExplorerTest() {
//
// 	}
// };
//
// TEST_F(NodeExplorerTest, TestConstruction) {
// 	ASSERT_TRUE(explorer.feasibilityCuts.empty());
// }
//
// class CutTest: public testing::Test {
// 	// using namespace Inavap;
//
// protected:
// 	// using namespace Inavap;
//
// 	// Inavap::Cut cut;
// 	//
// 	// vector<pair<uint64_t, double>> coeff = {
// 	// 	{Inavap::getKey(2,3,5), 43}, {Inavap::getKey(2,3,6), 90},
// 	// 	{Inavap::getKey(2,4,5), 5345}
// 	// };
//
// 	Inavap::Cut cut1{12, {
// 			{Inavap::getKey(2,3,5), 43}, {Inavap::getKey(2,3,6), 90},
// 		{Inavap::getKey(2,4,5), 5345}, {getKey(0,0,0), -4324}
// 	} };
//
// 	Inavap::Cut cut {-9821, { {281539401220100,306},
// 		{281603825729543,2664},{562997198454798,2790},{42950066190,2328},{281513632399378,663},
// 		{281513632530450,663},{281513632006162,663},{563061624012827,422},{171800133659,1580},
// 		{281509338284057,462},{281509336645657,231},{563074509373456,450},{51541508112,173},{563074507931664,450},
// 		{51540066320,173},{563074509504528,450},{51541639184,173},{1125981511352333,506},{38654836749,426},
// 		{64424640525,192},{47244771341,3132},{844506536214541,314},{38656409613,234},{47246344205,2940},
// 		{844506535297037,314},{38655492109,234},{47245426701,2940} }};
//
// 	CutTest() {
// 		cut.printOffsets();
// 	}
// };
//
// TEST_F(CutTest, TestCutKey) {
// 	uint64_t key = getKey(4,0,15);
// 	double actual = cut.get(key);
// 	ASSERT_EQ(actual, 306);
// 	key = getKey(7,0,30);
// 	actual = cut.get(key);
// 	ASSERT_EQ(actual, 2664);
// 	key = getKey(13,2,9);
// 	actual = cut.get(key);
// 	ASSERT_EQ(actual, 426);
// 	key = getKey(16,31,29);
// 	actual = cut.get(key);
// 	ASSERT_EQ(actual, 450);
// }
//
// TEST_F(CutTest, TestCutGetKey) {
//
// 	auto key = Inavap::getKey(2,4,55);
// 	double val = cut.get(key);
// 	cout << "Val from get: " << val << endl;
// 	ASSERT_EQ(val, 0);
// }
//
// TEST_F(CutTest, TestGetStart) {
// 	// cut.printOffsets();
// 	uint32_t key = Inavap::getKey(0,0,6) & 0XFFFFFFFF;
// 	ASSERT_EQ(cut.getStart(key), 2);
// }
//
// TEST_F(CutTest, TestWrongKey) {
// 	{
// 		auto key = getKey(14,6,10);
// 		ASSERT_EQ(cut.get(key), 2328);
// 		key = getKey(655,342,2);
// 		ASSERT_EQ(cut.get(key), 0);
// 		key = getKey(13,13,4);
// 		ASSERT_EQ(cut.get(key), 0);
// 		key = getKey(25,28,8);
// 		ASSERT_EQ(cut.get(key), 462);
// 		key = getKey(18,17,9);
// 		ASSERT_EQ(cut.get(key), 663);
// 		key = getKey(0,0,0);
// 		ASSERT_EQ(cut.get(key), 0);
// 	}
// 	// {
// 	// 	auto key = getKey(4,56,4);
// 	// 	double val = cut.get(key);
// 	// 	// cout << "Val from get: " << val << endl;
// 	// 	ASSERT_EQ(val, 0);
// 	// }
// 	//
// 	// {
// 	// 	auto key = getKey(2,3,6);
// 	// 	double val = cut.get(key);
// 	// 	// cout << "Val from get: " << val << endl;
// 	// 	ASSERT_EQ(val, 90);
// 	// }
// 	// {
// 	// 	auto key = getKey(0,0,0);
// 	// 	double val = cut.get(key);
// 	// 	ASSERT_EQ(val, -4324);
// 	// }
// }
//
// class DDTest: public testing::Test {
// protected:
// 	const string fileName {"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_91_20_1.txt"};
// 	const shared_ptr<Network> networkPtr = make_shared<Network>(Network{fileName});
// 	DD relaxedDD {networkPtr,EXACT};
// 	DD restrictedDD {networkPtr,RESTRICTED};
// 	Inavap::RestrictedDD new_restrictedDD {networkPtr,MAX_WIDTH};
// 	Inavap::RelaxedDD new_relaxedDD{networkPtr};
// 	vector<Inavap::Node> new_exact_cutset;
// 	vector<Node_t> old_exact_cutset;
//
// 	vector<pair<uint64_t, double>> coeff = {
// 		/*{Inavap::getKey(2,3,6), 90}, */
// 		{Inavap::getKey(1,2,3), 40.0}, {Inavap::getKey(1,2,4), 3.0}
// 	};
// 	Inavap::Cut optimalityCut{23.45, coeff};
//
// 	DDTest() {
// 		DDNode root{0};
// 		root.nodeLayer = 0;
// 		root.globalLayer = 0;
// 		string message =  "Compiled Relaxed Tree in ";
// 		// MEASURE_EXECUTION_TIME(relaxedDD.build(root),message);
// 		message = "Compiled Restricted Tree in ";
// 		old_exact_cutset = restrictedDD.build(root).value();
// 		cout << "Compiled old restricted Tree" << endl;
// 		relaxedDD.build(root);
// 		// MEASURE_EXECUTION_TIME(restrictedDD.build(root), message);
// 		//cout << "Network processing size: " << network.processingOrder.size() << endl;
// 		//cout << "Relaxed DD size: " << relaxedDD.tree.size() << endl;
// 		//cout << "Restricted DD size: " << restrictedDD.tree.size() << endl;
//
// 		Inavap::Node node{};
// 		Inavap::Node node2{};
// 		cout << "Compiling Restricted Tree. " << endl;
// 		new_exact_cutset = new_restrictedDD.buildTree(Inavap::Node()).value();
// 		cout << "Compiling relaxed tree " << endl;
// 		new_relaxedDD.buildTree(Inavap::Node());
// 	}
// };
//
// TEST_F(DDTest, TestDDCutsetSize) {
// 	// Size of cutset should be same. Solution vectors might change.
// 	ASSERT_EQ(new_exact_cutset.size(), old_exact_cutset.size());
// 	// compare the sizes of solution vector of each node.
// 	for (int i = 0; i < new_exact_cutset.size(); i++) {
// 		ASSERT_EQ(new_exact_cutset[i].solutionVector.size(), old_exact_cutset[i].solutionVector.size());
// 	}
// }
//
// TEST_F(DDTest, TestTreeSizes) {
// 	ASSERT_EQ(restrictedDD.tree.size(), new_restrictedDD.tree.size());
// 	ASSERT_EQ(relaxedDD.tree.size(), new_relaxedDD.tree.size());
// 	ASSERT_EQ(restrictedDD.tree.size(), relaxedDD.tree.size());
// 	// test layer sizes of all trees.
// 	for (int i = 0; i < restrictedDD.tree.size(); i++) {
// 		ASSERT_EQ(restrictedDD.tree[i].size(), new_restrictedDD.tree[i].size());
// 		ASSERT_EQ(relaxedDD.tree[i].size(), new_relaxedDD.tree[i].size());
// 	}
// }
//
// TEST_F(DDTest, TestRelaxedDD) {
// 	// check if layer sizes in corresponding trees are same or not.
// 	cout << "Tree: Relaxed DD" << endl;
// 	for (int i = 0; i < relaxedDD.tree.size(); i++) {
// 		cout << relaxedDD.tree[i].size() << " ";
// 	}cout << endl;
// 	cout << "Tree : new relaxed dd: "<< endl;
// 	for (int i = 0; i < relaxedDD.tree.size(); i++) {
// 		cout << new_relaxedDD.tree[i].size() << " ";
// 	}cout << endl;
// 	for (int i = 0; i < restrictedDD.tree.size(); i++) {
// 		ASSERT_EQ(restrictedDD.tree[i].size(), new_restrictedDD.tree[i].size());
// 	}
// }
//
// TEST_F(DDTest, TestDDRestrictedOptimalityCut) {
// 	double val = new_restrictedDD.applyOptimalityCut(optimalityCut);
// 	cout << "Optimality cut value is " << val << endl;
// }
//
// TEST_F(DDTest, TestRelaxedD) {
// 	// relaxed.
// }


// test with some constant data.
TEST(ConstantCut, SHOULD_RETURN_SAME) {
	map<tuple<int,int,int>, double> coeff;
	coeff.insert(make_pair(make_tuple(2,30,123), 432.0));
	coeff.insert(make_pair(make_tuple(2,30,124), 456.67));
	coeff.insert(make_pair(make_tuple(1,18,123), 1234.56));
	coeff.insert(make_pair(make_tuple(4,1,90), -1298.98));
	coeff.insert(make_pair(make_tuple(5,6,7), -1298.98));

	::Cut cut {FEASIBILITY, 3012.0321, coeff};

	Inavap::Cut cut2 = cutToCut(cut, nullptr);

	uint64_t key = getKey(static_cast<uint64_t>(30), static_cast<uint64_t>(2), static_cast<uint64_t>(123));
	ASSERT_EQ(cut2.get(key), 432.0);
	ASSERT_EQ(cut2.get(key), cut.get(2,30,123));
	key = getKey(static_cast<uint64_t>(30), static_cast<uint64_t>(2), static_cast<uint64_t>(124));
	ASSERT_EQ(cut2.get(key), 456.67);
	key = getKey(static_cast<uint64_t>(6), static_cast<uint64_t>(5), static_cast<uint64_t>(7));
	ASSERT_EQ(cut2.get(key), -1298.98);
	key = getKey(static_cast<uint64_t>(6), static_cast<uint64_t>(5), static_cast<uint64_t>(1));
	ASSERT_EQ(cut2.get(key), 0.0);

}

::Cut returnCut() {
	// generate artifical cut with
	map<tuple<int,int,int>,double> cut;

	std::random_device r;
	std::default_random_engine e(r());
	std::uniform_real_distribution<double> urd(0,100);
	std::uniform_int_distribution<int> uid(0,10);

	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < 10; j++) {
			for (int k = 0; k < 10; k++) {
				if (int num = uid(e); num %2 == 0) {
					auto val = urd(e);
					cut[make_tuple(i,j,k)] = val;
					// cout << val << endl;
				}
				else cut[make_tuple(i,j,j)] = 0;
			}
		}
	}

	double rhs = urd(e)*5;
	return ::Cut{FEASIBILITY, rhs, cut};
}

::Cut getActualCut(const Network* network, CutType type= OPTIMALITY) { // builds pseudo-cut from the given network. builds at random.
	map<tuple<int,int,int>, double> coeffs;

	std::random_device r;
	std::default_random_engine e(r());
	std::uniform_real_distribution<double> urd(-100,100);
	std::uniform_int_distribution<int> uid(0,10);

	for (const auto [a,b] : network->processingOrder) {
		// b is arc id
		const NetworkArc& arc = network->networkArcs[b];
		uint i = arc.tailId;
		uint q = arc.headId;

		// j is the possible decisions.
		const NetworkNode& node = network->networkNodes[q];
		for (uint j : node.outNodeIds) {
			// random coefficient.
			if (int num = uid(e); num % 2 == 0)
				coeffs[tuple(i,q,j)] = urd(e);
			else coeffs[tuple(i,q,j)] = 0.0;
		}
	}
	double RHS = urd(e);
	if (type == OPTIMALITY) {RHS *= 10;}
	else RHS = abs(RHS)*7;
	return ::Cut{type, RHS, coeffs};
}

Inavap::Cut changeCut(::Cut& cut) {
	vector<pair<uint64_t, double>> coeff;
	for (const auto [k,v] : cut.cutCoeff) {
		if (v == 0) continue;
		// auto [i,q,j] = k;
		uint64_t i = std::get<0>(k);
		uint64_t q = std::get<1>(k);
		uint64_t j = std::get<2>(k);
		uint64_t key = getKey(q,i,j);
		coeff.push_back(make_pair(key,v));
	}
	std::reverse(coeff.begin(), coeff.end());
	return Inavap::Cut{cut.RHS, coeff};
}

class TestInavapCut: public testing::Test {
protected:
	::Cut oldCut;
	Inavap::Cut newCut;

	TestInavapCut() : oldCut{returnCut()}, newCut(changeCut(oldCut))  {
		// build InavapCut;
	}
};
TEST_F(TestInavapCut, INVARIANT_0) {
	// size of Inavap::cut should be atmost ::Cut
	ASSERT_LE(newCut.coeff.size(), oldCut.cutCoeff.size());
	ASSERT_EQ(newCut.getRHS(), oldCut.RHS); // RHS should match.
}
TEST_F(TestInavapCut, INVARIANT_1) {
	// number of non-zero coeffs match with size of inavap cut.
	size_t oldCutSize = oldCut.cutCoeff.size();
	size_t newCutSize = newCut.coeff.size();
	size_t count = 0;
	for (const auto [k,v] : oldCut.cutCoeff) {
		if (v != 0) count++;
	}
	ASSERT_EQ(count, newCutSize);
	ASSERT_EQ(oldCutSize, newCutSize + (oldCutSize-count));
}
TEST_F(TestInavapCut, INVARIANT_2) {
	// coeff should match
	for (const auto [k,v] : oldCut.cutCoeff) {
		// auto [i,q,j] = k;
		uint64_t i = std::get<0>(k);
		uint64_t q = std::get<1>(k);
		uint64_t j = std::get<2>(k);
		uint64_t key = getKey(q,i,j);
		double val = newCut.get(key);
		ASSERT_EQ(val, v);
	}
}


class TestRestrictedDD: public testing::Test {
protected:
	const string fileName {"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_91_20_1.txt"};
	const shared_ptr<Network> networkPtr = make_shared<Network>(Network{fileName});
	Inavap::RestrictedDDNew new_restrictedDD;
	Inavap::RestrictedDD restrictedDD;
	DD old_restrictedDD{networkPtr,RESTRICTED};
	vector<Inavap::Node> cutset;
	vector<Inavap::Node> new_cutset;
	vector<Node_t> old_cutset;

	TestRestrictedDD(): new_restrictedDD{networkPtr,32}, restrictedDD{networkPtr,32} {
		cout << "File: " << fileName << endl;
		Node root;
		// root.globalLayer = 2;
		// root.solutionVector = {-1,-1};
		// root.states = {-1};
		DDNode node;
		// node.globalLayer = 2;
		// node.solutionVector = {-1,-1};
		// node.states = {-1};
		new_cutset = new_restrictedDD.buildTree(root).value();
		// cutset = restrictedDD.buildTree(root).value();
		old_cutset = old_restrictedDD.build(node).value();
	}
};

/*
 * INVARIANT 0: number of nodes == terminalId+1;
 * INVARIANT 1 : no arc with id '0'.
 * INVARIANT 2 : one parent for all nodes except root and terminal nodes.
 * INVARIANT 3 : node.incomingArcId == node.id for all nodes except root and terminal.
 * INVARIANT 4 : number of terminal arcs == number of nodes in last layer.
 */
TEST_F(TestRestrictedDD, INVARIANT_0) {
	ASSERT_TRUE(new_restrictedDD.nodes.size() == new_restrictedDD.terminalId+1);
	ASSERT_TRUE(new_restrictedDD.nodes.size() > 2);
	cout << "Number of nodes in the restricted DD: " << new_restrictedDD.nodes.size() << endl;
	ASSERT_EQ(new_restrictedDD.nodes.size(),old_restrictedDD.nodes.size());
	ASSERT_EQ(new_restrictedDD.arcs.size(), old_restrictedDD.arcs.size());
}
TEST_F(TestRestrictedDD, INVARIANT_1) {
	ASSERT_FALSE(new_restrictedDD.arcs.contains(0));
	ASSERT_FALSE(old_restrictedDD.arcs.contains(0));
}
TEST_F(TestRestrictedDD, INVARIANT_2) {
	for (size_t i = 1; i < new_restrictedDD.terminalId; i++) {
		ASSERT_TRUE(new_restrictedDD.nodes[i].incomingArc != 0);
	}
}
TEST_F(TestRestrictedDD, INVARIANT_3) {
	for (size_t i = 1; i < new_restrictedDD.terminalId; i++) {
		const Inavap::RestrictedDDNew::RDDNode& node = new_restrictedDD.nodes[i];
		ASSERT_EQ(node.incomingArc, node.id);
		const ::DDNode& node2 = old_restrictedDD.nodes[i];
		ASSERT_EQ(node2.incomingArcs[0], node2.id);
	}
}
TEST_F(TestRestrictedDD, INVARIANT_4) {
	size_t tree_size = new_restrictedDD.tree.size();
	ASSERT_EQ(new_restrictedDD.tree[tree_size-2].size(), new_restrictedDD.terminalInArcs.size());
}
TEST_F(TestRestrictedDD, INVARIANT_5) { // cutset should contain nodes if tree is not exact.
	ASSERT_EQ(!new_cutset.empty(), !new_restrictedDD.isTreeExact());
	ASSERT_EQ(new_cutset.size(), old_cutset.size());

	// verify both cutsets contains same nodes.
	for (size_t i = 0; i < new_cutset.size(); i++) {
		// compare solutions of nodes in cutset.
		const Inavap::Node& node1 = new_cutset[i];
		const ::Node_t& node2 = old_cutset[i];

		ASSERT_EQ(node1.solutionVector.size(), node2.solutionVector.size());
		ASSERT_EQ(node1.globalLayer, node2.globalLayer);
		ASSERT_EQ(node1.states.size(), node2.states.size());
		for (size_t j = 0; j < node1.states.size(); j++) {
			ASSERT_EQ(node1.states[j], node2.states[j]);
		}
		for (size_t j = 0; j < node1.solutionVector.size(); j++) {
			ASSERT_EQ(node1.solutionVector[j], node2.solutionVector[j]);
		}
	}
}
TEST_F(TestRestrictedDD, INVARIANT_6) {
	// incoming arc's decision should not be in the node's states (unless state updated)
	// should go layer by layer.
	const unordered_map<uint, Inavap::RestrictedDDNew::RDDNode>& nodes = new_restrictedDD.nodes;
	const unordered_map<uint, Inavap::DDArc>& arcs = new_restrictedDD.arcs;
	for (size_t layer = 1; layer < new_restrictedDD.tree.size()-1; layer++) {
		if (networkPtr->hasStateChanged[layer]) continue;
		for (auto id: new_restrictedDD.tree[layer]) {
			const Inavap::RestrictedDDNew::RDDNode& node = nodes.at(id);
			const Inavap::DDArc& arc = arcs.at(node.incomingArc);
			auto res = find(node.states.begin(), node.states.end(), arc.decision);
			if (arc.decision != -1)
				ASSERT_EQ(res, node.states.end());
			else ASSERT_NE(res, node.states.end());
		}
	}
}

TEST_F(TestRestrictedDD, INVARIANT_7) {
	// no empty layers.
	ASSERT_TRUE(!new_restrictedDD.tree.empty());
	ASSERT_TRUE(!old_restrictedDD.tree.empty());
	for (const Layer& layer: new_restrictedDD.tree) {
		ASSERT_TRUE(!layer.empty());
	}
	for (size_t i = 0; i < new_restrictedDD.tree.size()-1; i++) {
		ASSERT_EQ(new_restrictedDD.tree[i].size(), old_restrictedDD.tree[i].size());
		for (size_t j = 0; j < new_restrictedDD.tree[i].size(); j++) {
			ASSERT_EQ(new_restrictedDD.tree[i][j], old_restrictedDD.tree[i][j]);
		}
	}
}

/* Okay the suspect must be optimality code then?. Build restricted tree exactly as ::DD(restricted) and
 * make sure nodes and arcs match in both trees. Apply cut on both trees,*/
void printOldCut(const ::Cut& cut) {

	for (const auto[k,v] : cut.cutCoeff) {
		if (v != 0) {
			auto [i,q,j] = k;
			cout << "[("<< i << "," <<q << "," << j << "):"<<v<<"]" << endl;
		}
	}
}

TEST_F(TestRestrictedDD, INVARIANT_8) {
	// apply 1000 optimality cuts.
	for (size_t i = 0; i < 10000; i++) {
		::Cut old_cut = getActualCut(networkPtr.get());
		// apply this cut
		// printOldCut(old_cut);
		Inavap::Cut new_cut = cutToCut(old_cut, networkPtr.get());
		// new_cut.printCut(0);
		double lowerBound = new_restrictedDD.applyOptimalityCut(new_cut);
		// cout << "Lower Bound from new DD: " << lowerBound << endl;
		double lb2 = old_restrictedDD.applyOptimalityCutRestrictedLatest(old_cut);
		// cout << "Lower Bound : " << lb2 <<" (old), "<< lowerBound <<" (new)" <<endl;
		ASSERT_EQ(lowerBound, lb2);
	}
}

TEST_F(TestRestrictedDD, INVARIANT_9) {
	// apply feasibility cuts.
	for (size_t i= 0; i < 1000; i++) {
		::Cut old_cut = getActualCut(networkPtr.get(), FEASIBILITY);
		Inavap::Cut new_cut = cutToCut(old_cut, networkPtr.get());
		uint8_t res = new_restrictedDD.applyFeasibilityCut(new_cut);
		uint8_t res2 = old_restrictedDD.applyFeasibilityCutRestrictedLatest(old_cut);
		if (res != res2) {
			cout << "Failed in " << i << " th iteration" << endl;
		}
		ASSERT_EQ(res, res2);
		if (!res) {
			cout << "Breaking the loop " << endl;
			break;
		}
		cout << "Passed " << i << " iterations" << endl;
	}
}

TEST_F(TestRestrictedDD, INVARIANT_10) {
	// apply random cuts
	srand(time(NULL));
	for (size_t i = 0; i< 10; i++) {
		if (int r = rand(); r % 2 == 0) {
			// feasibility cut
			::Cut old_cut = getActualCut(networkPtr.get(), FEASIBILITY);
			Inavap::Cut new_cut = cutToCut(old_cut, networkPtr.get());
			uint8_t res = new_restrictedDD.applyFeasibilityCut(new_cut);
			uint8_t res2 = old_restrictedDD.applyFeasibilityCutRestrictedLatest(old_cut);
			ASSERT_EQ(res, res2);
			if (res == res2 && res2 == 0) {cout << "Breaking loop: "<< endl; break;}
		}
		else {
			// optimality cut
			::Cut old_cut = getActualCut(networkPtr.get(), OPTIMALITY);
			Inavap::Cut new_cut = cutToCut(old_cut, networkPtr.get());
			double lb = new_restrictedDD.applyOptimalityCut(new_cut);
			double lb2 = old_restrictedDD.applyOptimalityCutRestrictedLatest(old_cut);
			ASSERT_EQ(lb, lb2);
		}
	}
}


class TestRelaxedDD : public ::testing::Test {
protected:
	const string fileName {"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_91_20_1.txt"};
	const shared_ptr<Network> networkPtr = make_shared<Network>(Network{fileName});
	Inavap::RelaxedDD relaxedDD;
	::DD old_relaxedDD;

	// don't need cutset.
	TestRelaxedDD(): relaxedDD{networkPtr}, old_relaxedDD{networkPtr, EXACT} {
		Node root;
		::DDNode rootNode;
		relaxedDD.buildTree(root);
		old_relaxedDD.build(rootNode);
	}

};

TEST_F(TestRelaxedDD, TEST_NUMBER_OF_ARCS_MATCH) {
	ASSERT_EQ(relaxedDD.nodes.size(), old_relaxedDD.nodes.size());
}

int main() {
	testing::InitGoogleTest();
	return RUN_ALL_TESTS();
}