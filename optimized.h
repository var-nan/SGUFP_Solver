//
// Created by nandgate on 12/3/2024.
//

#ifndef OPTIMIZED_H
#define OPTIMIZED_H

#include "Network.h"
#include "Cut.h"
#include <mutex>
#include <shared_mutex>

const size_t FEASIBILITY_CONTAINER_CAPACITY = 100;
const size_t OPTIMALITY_CONTAINER_CAPACITY = 100;

namespace Inavap {

    using namespace std;

    typedef shared_lock<shared_mutex> ReaderLock;
    typedef unique_lock<shared_mutex> WriterLock;

    class CutResource {
        std::atomic<uint> count = 0; // add padding to different cache line. // sum of opt and feas containers.
        shared_mutex m;
        // add padding to put in different cacheline.
        vector<CutContainer *> fCutContainers;
        vector<CutContainer *> oCutContainers;

    public:
        CutResource() = default;

        [[nodiscard]] uint getCount() const noexcept {return count.load(memory_order::relaxed);}
        void add(pair<vector<CutContainer *>, vector<CutContainer *>> cuts);
        [[nodiscard]] pair<vector<CutContainer *>, vector<CutContainer *>> get(uint nF, uint nO);

        CutResource(const CutResource &other) = delete;
        CutResource(CutResource &&other) = delete;
        CutResource& operator=(const CutResource &other) = delete;
        CutResource& operator=(CutResource &&other) = delete;
    };
}




#endif //OPTIMIZED_H
