//
//  Common.h
//  dbfe
//
//  Created by Nigel Delaney on 7/10/16.
//
//

#ifndef Common_h
#define Common_h

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <vector>
#include <limits>

// Calls Rf_runif
// see Rmath.h in RCpp source
#define RUNIF R::runif(0, 1)
#define RPOIS R::rpois

typedef std::vector<double> dvector;



// Forward declaration of classes to avoid problems
namespace dbfe {
  class DiscretizedDFE;
  class ObservedWell;
  class PopulationSize;
  class MutationCounter;
  struct TimeFitnessClass
  {
  public:
    double time;
    int Class;
  };
}

#endif /* Common_h */
