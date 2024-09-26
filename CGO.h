//
// Created by erfank on 9/24/2024.
//

#ifndef CGO_H
#define CGO_H
#include <vector>
#include <map>
#include "Network.h"

using namespace std;
class Cut {
    public:
    map<tuple<int, int, int>, double> cutCoef;
    double RHS;
    char type = 'N';
    // Cut(vector<vector<vector<float>>> &y , float rhs , char t ) {
    //     cutCoef;
    //     RHS=rhs;
    //     type=t;
    // }
};
void generateCut(vector<int> &W_solution, Network &network ,Cut &newCut, int s) ;


#endif //CGO_H
