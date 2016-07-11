#include <memory>

#include "Common.h"
#include "PopulationSize.h"
#include "MutationCounter.h"

namespace dbfe {
/* This is an observed unit of data.  */
class ObservedWell {
  public:
    std::shared_ptr<PopulationSize> PopSize;
    double TotalTransfers;
    double ObservedFitness;
    double AmountOfTimeLastRun;
    int NumberOfSimulationsLastRun = 0;
    /// <summary>
    /// A class that keeps track of all the mutations that have been observed in this well from each bin.
    /// Used to update the prior.
    /// </summary>
    MutationCounter MutCounter;
    /// <summary>
    /// 0 is the neutral class.
    /// </summary>
    const int BinClass;
    double TotalGenerations () const
    {
      return TotalTransfers * PopSize->GenerationsInBetweenTransfers;
    }
  ObservedWell(double Transfers, double W, const DiscretizedDFE& dfe, std::shared_ptr<PopulationSize> size);
};
}; // close namespace
