//
// Created by nandgate on 10/21/2024.
//

#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "../DD.h"
#include <chrono>
#include <gtest/gtest.h>

using namespace std;

#define MEASURE_EXECUTION_TIME(block,message) \
	{		\
		auto start = std::chrono::high_resolution_clock::now();		\
		block;														\
		auto end = std::chrono::high_resolution_clock::now();		\
		std::chrono::duration<double> duration = end - start;	\
		cout << message << duration.count() << " seconds" << endl;	\
	}


inline size_t getNumNodesDD(const DD& dd) {
	size_t count = 0;

	for (const auto& layer: dd.tree) {
		count += layer.size();
	}
	return count;
}

inline void printTreeStatistics(const DD& dd) {
	cout << "Number of layers in tree: " << dd.tree.size() << endl;
}

inline void printTree(const DD& dd) {
	for (const auto& layer: dd.tree) {
		for (const auto id: layer) cout << id << " ";
		cout << endl;
	}
}
#endif //TEST_UTILS_H
