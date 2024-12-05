//
// Created by nandgate on 12/3/2024.
//

#include "optimized.h"

void Inavap::CutResource::add(pair<vector<CutContainer *>, vector<CutContainer *>> cuts) {
   // get pointer somehow? store them in heap. only master can update the structures.
   {
      WriterLock l{m};
      if (!cuts.first.empty()) fCutContainers.insert(fCutContainers.end(), cuts.first.begin(), cuts.first.end());
      if (!cuts.second.empty()) oCutContainers.insert(oCutContainers.end(), cuts.second.begin(), cuts.second.end());
   }
   count.store(fCutContainers.size() + oCutContainers.size(), memory_order::release);
}

/**
 * Returns the pair of feasibility cut containers and optimality cut containers from the global space.
 * The calling worker should maintain the invariant (nF+nO) <= count.
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