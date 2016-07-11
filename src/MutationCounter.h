#ifndef MutationCounter_h
#define MutationCounter_h

#include "Common.h"

namespace dbfe {

class MutationCounter {
public:
  std::vector<int> CountOfEachMutation;
  MutationCounter(const DiscretizedDFE& dfe);
  void AddMutations(std::vector<TimeFitnessClass>& Mutations);
  void AddMutationCounter(MutationCounter& MC);
  void AddCountToClass(int numberOfMutants, int MutantClass);
};
}; // namespace dbfe

#endif /* MutationCounter_h */
