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
        vector<Inavap::CutContainer *> fCutContainers;
        vector<Inavap::CutContainer *> oCutContainers;

    public:
        CutResource() = default;

        /** Returns the sum of sizes of feasibility cut containers and optimality cut containers. The memory order
         * of the load is set to relaxed.
         * <b>The caller should not derive any ordering with this function call.</b>
         */
        [[nodiscard]] uint getCount() const noexcept {return count.load(memory_order::relaxed);}
        void add(pair<vector<Inavap::CutContainer *>, vector<Inavap::CutContainer *>> cuts);
        [[nodiscard]] pair<vector<Inavap::CutContainer *>, vector<Inavap::CutContainer *>> get(uint nF, uint nO);

        CutResource(const CutResource &other) = delete;
        CutResource(CutResource &&other) = delete;
        CutResource& operator=(const CutResource &other) = delete;
        CutResource& operator=(CutResource &&other) = delete;
    };
}




#endif //OPTIMIZED_H
