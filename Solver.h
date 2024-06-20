//
// Created by nandgate on 6/3/24.
//
#pragma once

#include "Network.h"
#include "DD.h"
#include "grb.h"
#include <queue>

class Solver{

public:
	Network network;

	Solver(Network&& network1): network{network1}{};


	void solveStochasticDD();


};
