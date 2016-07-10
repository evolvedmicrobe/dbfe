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
#include <stdlib.h>
#include <string>
#include <vector>
#include <limits>


// Calls Rf_runif
// see Rmath.h in RCpp source
#define RUNIF R::runif(0, 1);
#define RPOIS R::rpois

typedef std::vector<double> dvector;

#endif /* Common_h */
