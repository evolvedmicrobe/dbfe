#include "Common.h"
namespace dbfe {
class PopulationSize {
public:
  double N0;
  double NF;
  double GenerationsInBetweenTransfers;
  double TotalGrowth() const {
    return (NF - N0);
  };
  double Ne() const {
    return N0 * log(2.0) * GenerationsInBetweenTransfers;
  };

  PopulationSize(double n0, double nf) : N0(n0), NF(nf)
  {
      GenerationsInBetweenTransfers = log2(nf / n0);
  };
};
} // close namespace
