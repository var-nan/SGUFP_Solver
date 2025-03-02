//
// Created by nandgate on 12/3/2024.
//

#ifndef OPTIMIZED_H
#define OPTIMIZED_H

#include "Network.h"
#include "Cut.h"
#include <mutex>
#include <atomic>
#include <shared_mutex>

#ifndef FEASIBILITY_CONTAINER_CAPACITY
    #define FEASIBILITY_CONTAINER_CAPACITY 12
#endif

#ifndef OPTIMALITY_CONTAINER_CAPACITY
    #define OPTIMALITY_CONTAINER_CAPACITY 12
#endif

namespace Inavap {

    using namespace std;

    typedef shared_lock<shared_mutex> ReaderLock;
    typedef unique_lock<shared_mutex> WriterLock;

    class CutResource {
        std::atomic<uint> count = std::atomic<uint>(0); // add padding to different cache line. // sum of opt and feas containers.
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

        void printStatistics() const noexcept {
            cout << "Number of Feasibility Cut containers " << fCutContainers.size() << endl;
            cout << "Number of Optimality Cut container " << oCutContainers.size() << endl;
        }

        /* TODO: create destructor that frees all the cut containers */
        ~CutResource() {
            // locking mutex not needed. Error if needed.
            // clear all the cuts from the memory.
            for (auto fCutContainer : fCutContainers) delete fCutContainer;
            for (auto oCutContainer : oCutContainers) delete oCutContainer;
        }
    };
}




#endif //OPTIMIZED_H
