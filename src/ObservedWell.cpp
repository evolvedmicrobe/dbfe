#include "ObservedWell.h"
#include "DiscretizedDFE.h"
#include "MutationCounter.h"
namespace dbfe {
ObservedWell::ObservedWell(double Transfers, double W, const DiscretizedDFE& dfe, std::shared_ptr<PopulationSize> size) :
  TotalTransfers(Transfers),
  ObservedFitness(W),
  PopSize(size),
  MutCounter(dfe),
  BinClass(dfe.AssignFitnessToBin(W))
{
  if (BinClass != 1)
  { MutCounter.AddCountToClass(1, BinClass); }
}
}; // close namespace
