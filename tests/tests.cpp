//
// Created by nandgate on 10/18/24.
//

#include <gtest/gtest.h>
#include "../DD.h"

using namespace std;

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
	const string fileName {"/home/nandgate/CLionProjects/SGUFP_Solver/40_50_1.txt"};
	Network network{fileName};
	DD relaxedDD {RELAXED};
	DD restrictedDD {RESTRICTED};

	DDTest() {
		DDNode root{0};
		relaxedDD.build(network, root, 0);
		restrictedDD.build(network, root, 0);

		cout << "Network processing size: " << network.processingOrder.size() << endl;
		cout << "Relaxed DD size: " << relaxedDD.tree.size() << endl;
		cout << "Restricted DD size: " << restrictedDD.tree.size() << endl;
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

//int main(){
//	testing::InitGoogleTest();
//	return RUN_ALL_TESTS();
//}