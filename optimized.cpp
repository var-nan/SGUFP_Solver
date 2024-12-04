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
}


pair<vector<CutContainer *>, vector<CutContainer *>> Inavap::CutResource::get(uint nF, uint nO){
   vector<CutContainer *> feas, opti;
   // reserve later.
   {
      ReaderLock l{m};
      if ( nF < fCutContainers.size()) {
         // copy
         for(; nF < fCutContainers.size(); nF++) feas.push_back(fCutContainers[nF]);
      }
      if (nO < oCutContainers.size()) {
         // copy pointers
         for (; nO < oCutContainers.size(); nO++) opti.push_back(oCutContainers[nO]);
      }
   }
   return {feas, opti};
}