#include "Common.h"

// A faster implementation of rmultinom accounting for the
// fact that I am only doing one draw, and am not validating data.
int multinomial(const dvector& probs) {
  auto rv = RUNIF;
  double cumSum = probs[0];
  int i = 0;
  do
  {
    if (rv < cumSum) { return i; }
    i++;
    cumSum += probs[i];

  }
  while (i < probs.size());
  Rcpp::stop("Multinomial sampling didn't happen correctly!");
}
