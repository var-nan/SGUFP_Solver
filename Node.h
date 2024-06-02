/**
 * Created by nandgate on 6/1/24
 */

#ifndef SGUFP_SOLVER_NODE_H
#define SGUFP_SOLVER_NODE_H

#include <vector>
#include <iostream>

typedef uint64_t State;

class Node {
public:
    uint32_t nodeId; /* id of the node */
    uint64_t objectiveVal; /* objective value is an unsigned long number */
    std::vector<uint32_t> state; /* node might have multiple state values when refining */
    uint32_t layerNo; /* with respect to global ordering of the variables */
    bool isExact;


    Node(uint32_t st, uint16_t layerNumber, uint64_t objective)
        : state{st}, layerNo{layerNumber}, objectiveVal{objective}, isExact{true} {}

    Node(const Node& node)
        : state{node.state}, layerNo{node.layerNo}, objectiveVal{node.objectiveVal}, isExact{node.isExact} {}

    /*
     * operators
     */

    Node& operator=(const Node& node) {

        return *this;
    }

    Node& operator=(Node&& node) noexcept {
        this->layerNo = node.layerNo;
        this->objectiveVal = node.objectiveVal;
        this->isExact = node.isExact;
        this->state = node.state;
        return *this;
    }

    bool operator==(const Node& node){
        return true; /* TODO: correct this ASAP */
    }

    bool operator<(const Node& node){
        return true; /* TODO: correct this ASAP */
    }

    bool operator<=(const Node& node){
        return true; /* TODO: correct this ASAP */
    }

    bool operator>(const Node& node){
        return true; /* TODO: correct this ASAP */
    }

    bool operator>=(const Node& node){
        return true; /* TODO: correct this ASAP */
    }

    struct HashFunction {
        size_t operator()(const Node& node) const {
            return std::hash<int>()(static_cast<int>(node.nodeId));
        }
    };
};

#endif //SGUFP_SOLVER_NODE_H
