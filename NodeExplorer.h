//
// Created by nandgate on 10/27/2024.
//

#ifndef NODEEXPLORER_H
#define NODEEXPLORER_H

#include <utility>

#include "DD.h"
#include "Cut.h"
#include "gurobi_c++.h"
#include "grb.h"
// #include "DDSolver.h"

// extern const Network network;

//

class OutObject{
public:
    double lb;
    double ub;
    vector<Node_t> nodes;
    bool success;

    OutObject(double lb_, double ub_, vector<Node_t> nodes_, bool s) :
        lb{lb_}, ub{ub_}, nodes{std::move(nodes_)}, success{s}{}
};

class NodeExplorer {

    //shared_ptr<Network> networkPtr;

    //CutContainer feasibilityCuts;
    //CutContainer optimalityCuts;
    GRBEnv env = GRBEnv();

    const shared_ptr<Network> networkPtr;

public:

    CutContainer feasibilityCuts;
    CutContainer optimalityCuts;

    explicit NodeExplorer(const shared_ptr<Network>& networkPtr_) : networkPtr{networkPtr_}, feasibilityCuts{FEASIBILITY}, optimalityCuts{OPTIMALITY} {
        // env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag,0);
        env.set(GRB_IntParam_Threads,1);
    }

    NodeExplorer(const shared_ptr<Network>& networkPtr_, pair<CutContainer, CutContainer> cuts): networkPtr{networkPtr_}, feasibilityCuts{cuts.first}, optimalityCuts{cuts.second} {
        env.set(GRB_IntParam_OutputFlag,0);
        env.set(GRB_IntParam_Threads,1);
    }

	pair<CutContainer, CutContainer> getCuts() noexcept { return {feasibilityCuts, optimalityCuts};}


    OutObject process(Node_t node, double optimalLB);
    OutObject process2(Node_t node, double optimalLB);
    OutObject process3(Node_t node, double optimalLB);
    OutObject process4(Node_t node, double optimalLB, const pair<vector<CutContainer *>, vector<CutContainer *>>& globalCuts);

    void clearCuts();

    #ifdef SOLVER_STATS
    void displayCutStats() const noexcept {
        cout << "Number of feasibility cuts in the container " << feasibilityCuts.cuts.size() << endl;
        cout << "Number of optimality cuts in the container: " << optimalityCuts.cuts.size() << endl;
    }
    #endif

};


namespace Inavap {

    /* Return Object from the Node Explorer to the DDSolver. */
    using OutObject = struct OutObj{
        double lb = std::numeric_limits<double>::lowest();
        /* reason: new (empty) nodes should not be processed (or sent to queue) if not initialized. */
        double ub = std::numeric_limits<double>::lowest();
        vector<Node> nodes;
        uint16_t status; // status bits after processing the node.
        /* this will hold the status of the node explorer operation on the given node.
         * 0X0 - fail
         * 0X1 - success
         * 0X2 - Success (with empty nodes returned).
         */

        // OutObj() = default;
        OutObj(const double lb_, const double ub_, vector<Node> nodes_, const uint16_t status_) :
                lb{lb_}, ub{ub_}, nodes{std::move(nodes_)}, status{status_}{}
        //

    };

    static auto INVALID_OBJECT = OutObject{std::numeric_limits<double>::lowest(),
                                                std::numeric_limits<double>::lowest(),
                                                    {}, false};

    class NodeExplorer {
        GRBEnv env = GRBEnv();
        const shared_ptr<Network> networkPtr;

    public:
        CutContainer feasibilityCuts; // local feasibility cuts.
        CutContainer optimalityCuts; // local optimality cuts.
        // vector<CutContainer *> globalFCuts;
        // vector<CutContainer *> globalOCuts;

        explicit NodeExplorer(const shared_ptr<Network>& networkPtr_): networkPtr{networkPtr_} {
            env.set(GRB_IntParam_OutputFlag,0);
            env.set(GRB_IntParam_Threads,1);
        }

        OutObject process(Node node, double optimalLB,
            const vector<CutContainer *> &globalFCuts, const vector<CutContainer *> &globalOCuts);

        OutObject process2(Node node, double optimalLB);

    };
}


#endif //NODEEXPLORER_H
