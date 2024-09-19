//
// Created by nandgate on 6/26/24.
//

#include <gtest/gtest.h>
#include "DD.h"
#include <iostream>


TEST(HelloTest, BasicAssertions){
	EXPECT_EQ(100,100);
	EXPECT_NE(1,2);
}

TEST(NetworkTest, TestNetworkRead){
	const string fileName{"/home/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};
	Network network(fileName);
	// some file name
	//ASSERT_EQ(network.Vbar.size(), 20);
	ASSERT_EQ(network.Vbar[0], 15);

}

class DDTest: public testing::Test {

protected:
	DDTest() {
//		diagram.arcs.insert({1, {1,0,1,0}});
//		diagram.arcs.insert({2,{2, 0,2,1}});
//		diagram.arcs.insert({3, {3,0,3,2}});
//		diagram.arcs.insert({4, {4,1,4,0}});
//		diagram.arcs.insert({5, {5,2,5,1}});
//		diagram.arcs.insert({6, {6,3,6,2}});
//
//		diagram.nodes.insert({1,{1}}); diagram.nodes[1].incomingArcs.push_back(1); diagram.nodes[1].outgoingArcs.push_back(4);
//		diagram.nodes.insert({2,{2}}); diagram.nodes[2].incomingArcs.push_back(2); diagram.nodes[2].outgoingArcs.push_back(5);
//		diagram.nodes.insert({3,{3}}); diagram.nodes[3].incomingArcs.push_back(3); diagram.nodes[3].outgoingArcs.push_back(6);
//		diagram.nodes.insert({4,{4}}); diagram.nodes[4].incomingArcs.push_back(4);
//		diagram.nodes.insert({5,{5}});
//		diagram.nodes.insert({6,{6}});
//		diagram.nodes.insert({0,{0}}); diagram.nodes[0].outgoingArcs.push_back(1); diagram.nodes[0].outgoingArcs.push_back(2);
//		diagram.nodes[0].outgoingArcs.push_back(3);
//
//		diagram.lastInserted = 6;
		const string fileName{"/home/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};
		Network network{fileName};
		DDNode node{0};
		diagram.build(network, node, 0);
	}

	DD diagram;
};

TEST_F(DDTest, TestDeleteArc){
	ASSERT_EQ(diagram.arcs.size(), 6);
	diagram.deleteArcById(2);
	ASSERT_EQ(diagram.arcs.size(), 5);
	ASSERT_EQ(diagram.nodes[2].incomingArcs.size(), 0);
}

TEST_F(DDTest, TestBuild){
	const string fileName {"/home/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};
	Network network {fileName};

	DD dd;
	DDNode node{0};

	dd.build(network, node, 0);

	ASSERT_EQ(dd.tree.size(), 5);
	ASSERT_EQ(dd.tree[2].size(), 4);
	ASSERT_EQ(dd.tree[3].size(), 6);
	ASSERT_EQ(dd.tree[4].size(), 8);
	ASSERT_EQ(5, dd.exactLayer);
}

TEST_F(DDTest, TestMerge){

	DDNode& node1 = diagram.nodes[2];
	DDNode& node2 = diagram.nodes[3];
	size_t node1_incoming = node1.incomingArcs.size();
	size_t node2_incoming = node2.incomingArcs.size();

	//ASSERT_TRUE(diagram.nodes[3].incomingArcs.size());
	diagram.mergeNodes(node1, node2); // after merging, incoming arcs should be updated.
	ASSERT_EQ(diagram.nodes[3].incomingArcs.size(), node2_incoming-1);
	ASSERT_EQ(diagram.nodes[2].incomingArcs.size(), node1_incoming+1);
	// TODO: test more assumptions.
}

TEST_F(DDTest, TestBuildNextLayer){
	// create two layers. current and next.
	vector<int> currentLayer = {3,4,5,6};
	vector<int> nextLayer;

	diagram.updateState(currentLayer, {-1,2,4});

	diagram.buildNextLayer(currentLayer, nextLayer, 2);

	ASSERT_TRUE(!nextLayer.empty());

	for (const auto id: nextLayer){
		cout << id << " ";
	}
	cout << endl;


}

TEST_F(DDTest, TestDuplicate){
	size_t s_inc_size = diagram.nodes[2].incomingArcs.size();
	size_t t_inc_size = diagram.nodes[3].incomingArcs.size();
	diagram.mergeNodes(diagram.nodes[2], diagram.nodes[3]);
	ASSERT_EQ(diagram.nodes[2].incomingArcs.size(), s_inc_size+t_inc_size );
	diagram.nodes.erase(3); // only in this function.
	//s_inc_size = diagram.nodes[2].incomingArcs.size();
	diagram.duplicateNode(2);
	ASSERT_EQ(diagram.nodes[2].incomingArcs.size(), 1);
}

//TEST(DDTreeTest, TestDeleteArcById){
//	DD diagram{};
//	diagram.arcs.insert({1, {1,0,1,0}});
//	diagram.arcs.insert({2,{2, 0,2,1}});
//	diagram.arcs.insert({3, {3,0,3,2}});
//	diagram.arcs.insert({4, {4,1,4,0}});
//	diagram.arcs.insert({5, {5,2,5,1}});
//	diagram.arcs.insert({6, {6,3,6,2}});
//
//	diagram.nodes.insert({1,{1}}); diagram.nodes[1].incomingArcs.push_back(1); diagram.nodes[1].outgoingArcs.push_back(4);
//	diagram.nodes.insert({2,{2}}); diagram.nodes[2].incomingArcs.push_back(2); diagram.nodes[2].outgoingArcs.push_back(5);
//	diagram.nodes.insert({3,{3}}); diagram.nodes[3].incomingArcs.push_back(3); diagram.nodes[3].outgoingArcs.push_back(6);
//	diagram.nodes.insert({4,{4}}); diagram.nodes[4].incomingArcs.push_back(4);
//	diagram.nodes.insert({5,{5}});
//	diagram.nodes.insert({6,{6}});
//	diagram.nodes.insert({0,{0}}); diagram.nodes[0].outgoingArcs.push_back(1); diagram.nodes[0].outgoingArcs.push_back(2);
//	diagram.nodes[0].outgoingArcs.push_back(3);
//
//	ASSERT_EQ(diagram.arcs.size(), 6);
//	diagram.deleteArcById(2);
//	//std::cout << diagram.arcs.size() << std::endl;
//	ASSERT_EQ(diagram.arcs.size(), 5);
//	ASSERT_EQ(diagram.nodes[2].incomingArcs.size(), 0);
//	//
//}

bool compareStates(const unordered_set<int>& s1, const unordered_set<int>& s2){
	return s1 == s2;
}


TEST(DDTreeTest, TestGetSolutionFunction){
	// create a tree as a linked list.
	DD dd;

	DDNode rootNode{0};
	//rootNode.solutionVector = {};
	rootNode.solutionVector = {4,5,6,7};
	dd.startTree = rootNode.solutionVector.size();
	DDNode node1{1};
	DDNode node2{2};
	DDNode node3{3};
	DDNode node4{4};
	dd.exactLayer = 5;
	DDArc arc1{1, 0,1, -1};
	DDArc arc2{2, 1,2, 8};
	DDArc arc3{3, 2,3, 6};
	DDArc arc4{4, 3,4, 5};

	node1.incomingArcs.push_back(arc1.id);
	node1.outgoingArcs.push_back(arc2.id);
	node2.incomingArcs.push_back(arc2.id);
	node2.outgoingArcs.push_back(arc3.id);
	node3.incomingArcs.push_back(arc3.id);
	node3.incomingArcs.push_back(arc4.id);
	node4.incomingArcs.push_back(arc4.id);

	dd.nodes.insert(make_pair(0, rootNode));
	dd.nodes.insert(make_pair(1, node1));
	dd.nodes.insert(make_pair(2, node2));
	dd.nodes.insert(make_pair(3, node3));
	dd.nodes.insert(make_pair(4, node4));

	dd.arcs.insert(make_pair(1, arc1));
	dd.arcs.insert(make_pair(2, arc2));
	dd.arcs.insert(make_pair(3, arc3));
	dd.arcs.insert(make_pair(4, arc4));

	vi actualSolutionVector{4,5,6,7,-1, 8,6,5};

	auto solutionVector = dd.getSolutionVector(4);
	for (int i = 0; i < solutionVector.size(); i++) ASSERT_EQ(actualSolutionVector[i], solutionVector[i]);
}

TEST_F(DDTest, TestGetCutsetFunction){

	const auto cutset = diagram.getExactCutset();

	for (auto node: cutset){
		for (const auto i: node.solutionVector){ cout << i << " "; }
		cout << endl;
	}
}

TEST(DDTreeTest, TestStateUpdateFunction){
	DD ddtree;

	ddtree.nodes.insert({1,{1}});
	ddtree.nodes.insert({2,{2}});
	ddtree.nodes.insert({3, {3}});
	vector<int> currentLayer = {1,2,3};
	unordered_set<int> states = {6,5,4};
	ddtree.updateState(currentLayer, states);

	//ASSERT_TRUE(compareStates({2,4,1}, {1,2,3}));

	for (const auto id: currentLayer){
		ASSERT_TRUE(compareStates(ddtree.nodes[id].states, states));
	}
}

TEST(DDTreeTest, TestPruneNodeFunction){
	// test prune node function, build network, and constexpr.
	//RestrictedDD restrictedDD {16};
}

//int main(){
//	testing::InitGoogleTest();
//	return RUN_ALL_TESTS();
//}