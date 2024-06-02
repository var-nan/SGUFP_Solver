//
// Created by nandgate on 6/1/24.
//

#ifndef SGUFP_SOLVER_DD_H
#define SGUFP_SOLVER_DD_H

#include "Node.h"
#include <unordered_set>

typedef std::unordered_set<Node, Node::HashFunction> Layer;

class RelaxedDD{

private:
	uint32_t upperBound;
	uint32_t maxWidth;
	bool isExact;
	Layer cutset { };
	
	/*
	 * Merge nodes based on a selected criteria.
	 */
	void mergeNodes(Layer& layer) {
		/* TODO; complete this function */
		/* debug this function */
	}
	
	void insertNode(Layer& currentLayer, Node&& node){
	
	}
	
public:


};

class DD {

private:
    Node root;
    std::vector<uint32_t> coefficients;
    std::vector<uint32_t> objectives;

    uint64_t objective;

    void insertNode(Layer& currentLayer, const Node& node){
        /* insert node to given layer */

        // if node is already present in the layer, just modify its attributes.
        auto got = currentLayer.find(node);

        if (got == currentLayer.end()) /* node not present in the layer */
            currentLayer.insert(node);
        else {
            /* node not present in the layer */
            auto temp = currentLayer.extract(node).value();

            /* TODO: update attributes of the node with the given node */
        }
    }


public:

     DD(std::vector<uint32_t>&& coeff, std::vector<uint32_t>&& obj, uint32_t state)
        : coefficients{ coeff}, objectives{obj}, root{state, 0,0}, objective{0} {}

    uint64_t getObjectiveVal() {
        return this->objective;
    }

	void build(uint16_t startLayer=0){


	}

};


#endif //SGUFP_SOLVER_DD_H
