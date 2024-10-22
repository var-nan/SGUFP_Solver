//
// Created by nandgate on 10/18/24.
//

#include "test_utils.h"

TEST(test1, test1) {
	ASSERT_EQ(1,1);
}

TEST(test2, testShuffle) {

	size_t n = 1000;
	size_t m = 20;
	auto res = getShuffledList(n,m);
	ASSERT_EQ(res.size(), m);
	for (int i = 0; i < m-1; i++){
		ASSERT_LE(res[i], res[i+1]);
	}
}

class DDTest : public testing::Test {
protected:
	const string fileName {"C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};
	Network network{fileName};
	DD relaxedDD {RELAXED};
	DD restrictedDD {RESTRICTED};

	DDTest() {
		DDNode root{0};
		string message =  "Compiled Relaxed Tree in ";
		MEASURE_EXECUTION_TIME(relaxedDD.build(network, root, 0),message);
		message = "Compiled Restricted Tree in ";
		MEASURE_EXECUTION_TIME(restrictedDD.build(network, root, 0), message);
		//cout << "Network processing size: " << network.processingOrder.size() << endl;
		//cout << "Relaxed DD size: " << relaxedDD.tree.size() << endl;
		//cout << "Restricted DD size: " << restrictedDD.tree.size() << endl;
	}
};

TEST_F(DDTest, TestRestrictedDD) {
	for (const auto layer : restrictedDD.tree) {
		cout << layer.size() << " ";
	}
	cout << endl;
}

TEST_F(DDTest, TestRelaxedDD){
	for (const auto layer: relaxedDD.tree){
		cout << layer.size() << " ";
	}
	cout << endl;
}

TEST_F(DDTest, TestOutArcsAndState) {
	for (size_t i = 0; i < relaxedDD.tree.size()-2; i++){
		for (auto nodeId : relaxedDD.tree[i]) {
			const auto& node = relaxedDD.nodes[nodeId];
			ASSERT_EQ(node.outgoingArcs.size(), node.states.size());
		}
	}
}


TEST_F(DDTest, TestReduceLayer){
	for (int i = 1; i < restrictedDD.tree.size()-2; i++){
		int s1 = restrictedDD.tree[i].size();
		restrictedDD.reduceLayer(restrictedDD.tree[i]);
		int s2 = restrictedDD.tree[i].size();

		if (network.hasStateChanged[i])ASSERT_EQ(s2, 1);
		else ASSERT_GE(s1, s2);
	}
}


TEST_F(DDTest, TestDelete) {
	 // remove a node from the tree.
	auto id = restrictedDD.tree[1][0];
	auto count1 = getNumNodesDD(restrictedDD);
	cout << "Num nodes before delete " << count1 << endl;
	ASSERT_EQ(restrictedDD.nodes.size(), count1);
	auto start = std::chrono::high_resolution_clock::now();
	restrictedDD.removeNode(id);
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	auto count2 = getNumNodesDD(restrictedDD);
	cout << "Num nodes after delete " << count2 << endl;
	cout << "Removed " << (count1 - count2) << " nodes in " << elapsed_seconds.count() << " seconds." << endl;
	ASSERT_EQ(restrictedDD.nodes.size(), count2);
}

TEST_F(DDTest, TestBatchDelete) {
	auto id = restrictedDD.tree[1][0];
	auto id2 = restrictedDD.tree[1][1];
	auto count1 = getNumNodesDD(restrictedDD);
	ASSERT_EQ(restrictedDD.nodes.size(), count1);
	auto start = std::chrono::high_resolution_clock::now();
	restrictedDD.batchRemoveNodes({id,id2});
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	auto count2 = getNumNodesDD(restrictedDD);
	cout << "Num nodes after delete " << count2 << endl;
	cout << "Removed " << (count1 - count2) << " nodes in " << elapsed_seconds.count() << " seconds." << endl;
	ASSERT_EQ(restrictedDD.nodes.size(), count2);
}


TEST_F(DDTest, TestBatchDeleteFullLayer) {

	size_t totalNodes = getNumNodesDD(restrictedDD);
	// remove whole layer
	auto layer = restrictedDD.tree[2];
	ASSERT_TRUE(!restrictedDD.nodes.empty());
	const string message = "Removed " + to_string(totalNodes) + " nodes in ";
	MEASURE_EXECUTION_TIME(restrictedDD.batchRemoveNodes(layer), message);
	cout << " seconds." << endl;
	// tree should be empty.
	ASSERT_TRUE(restrictedDD.nodes.empty());
	ASSERT_TRUE(restrictedDD.arcs.empty());
}


class DDSubRootTest : public testing::Test {
protected:
	DD restrictedDD;
	const string fileName = "C:/Users/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt";
	Network network{fileName};
	DDSubRootTest() {


		DDNode root{3433};
		root.incomingArcs = {434,32};
		root.outgoingArcs = {3232,2323,323,23,434,12,43,890};
		root.states = {1,2,3,4,5,6,7,8,9}; // TODO: what happens to state in between?

		root.solutionVector = {4,3,6,7,87};
		const string message = "Compiled Restricted Tree in ";
		MEASURE_EXECUTION_TIME(restrictedDD.build(network, root, 2), message);

		printTreeStatistics(restrictedDD);
	}
};

TEST_F(DDSubRootTest, Test1) {
	ASSERT_TRUE(!restrictedDD.nodes[0].solutionVector.empty());
	ASSERT_TRUE(restrictedDD.nodes[0].incomingArcs.empty());
	// test number of layers is equal to number of processing arcs.
	auto ind = network.processingOrder.size()-restrictedDD.startTree;
	auto treeSize = restrictedDD.tree.size()-2;
	ASSERT_EQ(ind, treeSize);
	// test that node layers are properly assigned
	auto start = 0;
	for (auto& layer: restrictedDD.tree) {
		for (auto id: layer) {
			auto node = restrictedDD.nodes[id];
			ASSERT_EQ(node.nodeLayer, start);
		}
		start++;
	}
}

TEST_F(DDSubRootTest, Test2) {

	auto rootNode = restrictedDD.nodes[0];
	auto cutset = restrictedDD.getExactCutset();
	size_t pathSize = rootNode.solutionVector.size() + restrictedDD.exactLayer;
	for (const auto& node : cutset) {
		for (size_t i = 0; i < rootNode.solutionVector.size(); i++) // root's solution should be in node's solution.
			ASSERT_EQ(rootNode.solutionVector[i], node.solutionVector[i]);
		ASSERT_EQ(pathSize, node.solutionVector.size());
	}
}

// int main(){
// 	testing::InitGoogleTest();
// 	return RUN_ALL_TESTS();
// }