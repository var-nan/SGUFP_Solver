//
// Created by nandgate on 12/3/2024.
//

#include "optimized.h"

/**
 * Adds the new cuts to the global space. Only Master thread can call this function to
 * add new cuts to the global space after accumulating new cuts from all the Workers.
 *
 * Insertion of new cuts happen by appending the pointers (pointing to the cut container)
 * to the corresponding global cut container.
 * @param cuts - pair of vector of pointers to new cuts.
 */
void Inavap::CutResource::add(pair<vector<CutContainer *>, vector<CutContainer *>> cuts) {
   // get pointer somehow? store them in heap. only master can update the structures.
   {
      WriterLock l{m};
      if (!cuts.first.empty()) fCutContainers.insert(fCutContainers.end(), cuts.first.begin(), cuts.first.end());
      if (!cuts.second.empty()) oCutContainers.insert(oCutContainers.end(), cuts.second.begin(), cuts.second.end());
      count.store(fCutContainers.size() + oCutContainers.size(), memory_order::relaxed);
   }
   /* TODO: Some reader threads may observe large container than the count. Either handle this case from reader thread,
    * or move the store operation of count in the mutex block (might hit performance). */
   // assuming the optimizer/hardware scheduler DOESN'T change this order. TODO: change to relaxed store.
   // count.store(fCutContainers.size() + oCutContainers.size(), memory_order::release);
}

/**
 * Returns the pair of feasibility cut containers and optimality cut containers from the global space.
 * The calling worker should maintain the invariant (nF+nO) <= count.
 *
 * NOTE: It is possible that the function returns more cuts more than the count.
 * @param nF - number of feasibility cut containers the calling worker aware of.
 * @param nO - number of optimality cut containers the calling worker aware of.
 */
pair<vector<Inavap::CutContainer *>, vector<Inavap::CutContainer *>> Inavap::CutResource::get(uint nF, uint nO){
   vector<Inavap::CutContainer *> feas, opti;
   // reserve later.
   {
      ReaderLock l{m};
      while (nF < fCutContainers.size()) feas.push_back(fCutContainers[nF++]);
      while (nO < oCutContainers.size()) opti.push_back(oCutContainers[nO++]);
   }
   return {feas, opti};
}