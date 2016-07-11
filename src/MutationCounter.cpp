#include "MutationCounter.h"
#include "DiscretizedDFE.h"
namespace dbfe {

MutationCounter::MutationCounter(const DiscretizedDFE& dfe) :
    CountOfEachMutation(dfe.NumberOfClassesIncludingNeutral()) {}

void
MutationCounter::AddMutations(std::vector<TimeFitnessClass>& Mutations)
{
  for (auto tfc : Mutations)
  {
    CountOfEachMutation[tfc.Class] += 1;
  }
}
void
MutationCounter::AddCountToClass(int numberOfMutants, int MutantClass)
{
  CountOfEachMutation[MutantClass] += numberOfMutants;
}
void
MutationCounter::AddMutationCounter(MutationCounter& MC)
{

  if (MC.CountOfEachMutation.size() != CountOfEachMutation.size())
  { Rcpp::stop("Can't add mismatched mutant counters"); }
  //this.CountOfEachMutation = CountOfEachMutation.Zip(MC.CountOfEachMutation,(x,y)=>x+y).ToArray();
  for (int i=0; i < MC.CountOfEachMutation.size(); i++) {
    CountOfEachMutation[i] += MC.CountOfEachMutation[i];
  }
}

}
