//
// Created by nandgate on 6/1/24.
//

#ifndef SGUFP_SOLVER_DD_H
#define SGUFP_SOLVER_DD_H

#include <iostream>
#include <vector>

#include <unordered_set>

class DDNode {
public:
	uint32_t nodeId; /* id of the node */
	uint64_t objectiveVal; /* objective value is an unsigned long number */
	std::vector<uint32_t> state; /* node might have multiple state values when refining */
	uint32_t layerNo; /* with respect to global ordering of the variables */
	bool isExact;


	DDNode(uint32_t st, uint16_t layerNumber, uint64_t objective)
			: state{st}, layerNo{layerNumber}, objectiveVal{objective}, isExact{true} {}

	DDNode(const DDNode& node)
			: state{node.state}, layerNo{node.layerNo}, objectiveVal{node.objectiveVal}, isExact{node.isExact} {}

	/*
	 * operators
	 */

	DDNode& operator=(const DDNode& node) {

		return *this;
	}

	DDNode& operator=(DDNode&& node) noexcept {
		this->layerNo = node.layerNo;
		this->objectiveVal = node.objectiveVal;
		this->isExact = node.isExact;
		this->state = node.state;
		return *this;
	}

	bool operator==(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator<(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator<=(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator>(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	bool operator>=(const DDNode& node){
		return true; /* TODO: correct this ASAP */
	}

	struct HashFunction {
		size_t operator()(const DDNode& node) const {
			return std::hash<int>()(static_cast<int>(node.nodeId));
		}
	};
};


typedef std::unordered_set<DDNode, DDNode::HashFunction> Layer;

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
	
	void insertNode(Layer& currentLayer, DDNode&& node){
	
	}
	
public:


};

class DD {

private:
    DDNode root;
    std::vector<uint32_t> coefficients;
    std::vector<uint32_t> objectives;

    uint64_t objective;

    void insertNode(Layer& currentLayer, const DDNode& node){
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
