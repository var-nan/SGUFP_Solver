//
// Created by erfank on 9/27/2024.
//

#ifndef NEWCGO_H
#define NEWCGO_H
#include <iostream>
#include <vector>
#include <map>
#include "Network.h"

class Cut {
public:
    map<tuple<int, int, int>, double> cutCoef;
    double RHS;
    char type = 'N';
};

void newgenerateCut(vector<int> &W_solution, Network &network ,Cut &newCut, int s) ;



#endif //NEWCGO_H
