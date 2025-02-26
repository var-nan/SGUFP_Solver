//
// Created by nandgate on 12/3/2024.
//

#ifndef OPTIMIZED_H
#define OPTIMIZED_H

#include "Network.h"
#include "Cut.h"
// #include <mutex>
// #include <atomic>
// #include <shared_mutex>

#include <omp.h>

#ifndef FEASIBILITY_CONTAINER_CAPACITY
    #define FEASIBILITY_CONTAINER_CAPACITY 12
#endif

#ifndef OPTIMALITY_CONTAINER_CAPACITY
    #define OPTIMALITY_CONTAINER_CAPACITY 12
#endif

namespace Inavap {

    using namespace std;

    // typedef shared_lock<shared_mutex> ReaderLock;
    // typedef unique_lock<shared_mutex> WriterLock;

    class CutResource {

        /* cut resource is a vector of pointers of cut containers. */
        uint count = 0; // add padding to different cache line. // sum of opt and feas containers.
        // shared_mutex m;
        omp_lock_t lock;
        // add padding to put in different cacheline.
        vector<Inavap::CutContainer *> fCutContainers;
        vector<Inavap::CutContainer *> oCutContainers;

    public:
        CutResource() {
            omp_init_lock(&lock);
        }

        /** Returns the sum of sizes of feasibility cut containers and optimality cut containers.
         */
        [[nodiscard]] uint getCount() const noexcept {
            uint c  = 0;
            #pragma omp atomic read
            c = count;
            return c;
        }
        void add(pair<vector<Inavap::CutContainer *>, vector<Inavap::CutContainer *>> cuts);
        [[nodiscard]] pair<vector<Inavap::CutContainer *>, vector<Inavap::CutContainer *>> get(uint nF, uint nO);

        CutResource(const CutResource &other) = delete;
        CutResource(CutResource &&other) = delete;
        CutResource& operator=(const CutResource &other) = delete;
        CutResource& operator=(CutResource &&other) = delete;

        void printStatistics() const noexcept {
            cout << "Number of Feasibility Cut Containers " << fCutContainers.size() << endl;
            cout << "Number of Optimality Cut Containers " << oCutContainers.size() << endl;
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
