#include "Common.h"
#include "EvolvingPopulation.h"
#include <numeric>


namespace dbfe {


   inline dvector Element_A_Times_E_ToThe_RT(const dvector& A, const dvector&  r, const double t)
   {
        if (A.size() != r.size())
        { Rcpp::stop("Cannot Multiply vectors of Unequal Length"); }
        dvector toR(A.size());
        for (int i = 0; i < A.size(); i++)
        {
            toR[i] = A[i] * std::exp(r[i] * t) ;
        }
        return toR;
    }

    inline dvector ElementSubtract(const dvector& x, const dvector& y) {
    	dvector toR(x.size());
    	if (x.size() != y.size()) {
    		Rcpp::stop("Tried to subtract unequal size vectors");
    	}
    	for(int i=0; i < x.size(); i++) {
    		toR[i] = x[i] - y[i];
    	}
    }

	inline void ElementAddInPlace(dvector& x, const dvector& y)
    {
        if (x.size() != y.size() )
        { Rcpp::stop("Cannot Subtract and Sum Vectors of Unequal Length"); }
        for (int i = 0; i < x.size(); i++)
        {
            x[i] += y[i];
        }
    }

	EvolvingPopulation::EvolvingPopulation(spDiscretizedDFE dfe, const PopulationSize& popSize) :
	 dfe(dfe),
	 PopSizes(dfe->NumberOfClassesIncludingNeutral()),
	 popRelativeFitnesses(dfe->MidPoints),
	 MutCounter(*dfe)
	  {
          LogrelativePopIncrease = std::log(popSize.NF / popSize.N0);
          InvRelativePopIncrease = popSize.N0 / popSize.NF;
          PopSizes[0] = popSize.N0;
          N0 = popSize.N0;
          NF = popSize.NF;


	  }

void
EvolvingPopulation::GrowOneCycle() {

  /* CALCULATE MEAN FITNESS AND GROWTH TIMES */
        double curPopSize = std::accumulate(PopSizes.begin(), PopSizes.end(), 0.0);
        dvector freqs(PopSizes.size());
        std::transform(PopSizes.begin(), PopSizes.end(), freqs.begin(), [curPopSize](const double a) {return a / curPopSize;} );
        // Get average fitness
        double meanFitnessStart = std::inner_product(freqs.begin(), freqs.end(), popRelativeFitnesses.begin(), 0.0);
        //how much growth occurs, approximate by how long the mean fitness would take to grow
        double growthTime = std::log(NF / curPopSize) / meanFitnessStart;
        // Get the size of each population at the end
        dvector finalPopSizes = Element_A_Times_E_ToThe_RT(PopSizes, popRelativeFitnesses, growthTime);
        dvector NFminusN0akaTimeForPoisson  = ElementSubtract(finalPopSizes, PopSizes);//also the time for poisson.
        //Sum of N_f- N_0, note approximate as if a mutation occurs that one guy that mutates will no longer contribute to this time
        TotalTime += std::accumulate(NFminusN0akaTimeForPoisson.begin(), NFminusN0akaTimeForPoisson.end(), 0.0);


  /* CALCULATE HOW MANY MUTATIONS IN EACH CLASS */
        dvector newMuts(PopSizes.size());
        dvector& mutRates = dfe->MutationRates();
        for (int i = 0; i < PopSizes.size(); i++)
        {
            //Only work with populations that exist to mutate
            if (PopSizes[i] > 0)
            {
                //New algorithm, for each class we will go through and mutate
                // See: https://github.com/evolvedmicrobe/InferenceMachinary/commit/cdd03b289731bb28e687ee39941230e7755803a6
                //the expected number of each fitness class
                //TODO: The description in the write-up needs to be updated to account for this.
                //TODO: Decide later if in the case of very few classes (mean muts =0) we should not do one poisson for all classes and
                //assign after the fact. Note small poissons are ~3X more expensive than single unifoRm with the MersenneTwister RNG
                for(int j = i; j < mutRates.size(); j++) {
                  // Note population sizes includes the neutral class, while mutRates does not, so we start at the same value
                    //Determine the mean number of mutations
                    double meanMutForClass = mutRates[j] * NFminusN0akaTimeForPoisson[i];
                    //Sample poisson random
                    int mutNumber = RPOIS(meanMutForClass);
                    int w = j + 1; // Add one to convert from mutation class to fitness class
                    MutCounter.CountOfEachMutation[w] += mutNumber;
                    //If higher, add to population
                    if (mutNumber > 0)
                    {
                      // Figure out the max value on the time scale for uniform sampling
                      double maxRescaled = std::exp(growthTime * dfe->MidPoints[i]);
                      //Now generate a mutation for each class
                      for (int k = 0; k < mutNumber; k++)
                      {
                            double uniformTime = RUNIF * maxRescaled;
                            //Convert back to original scale.
                            double actualTime = std::log(uniformTime) / dfe->MidPoints[i]; //TODO: Verify that this is cleaned up by common sub-expression elimination
                            //Find out how much the mutant grew, then add it to
                            //the new group and subtract from the old population
                            //group the growth it would have contributed there
                            //had it not mutated
                            double gr = dfe->MidPoints[w];
                            double GT = std::exp(gr * (growthTime - actualTime));
                            newMuts[w] += GT;
                            double old_gr = dfe->MidPoints[i];
                            double oldGT = std::exp(old_gr * (growthTime - actualTime));
                            if (oldGT > finalPopSizes[i]) { //if the growth is higher due to the approximation, subtract it all, should rarely happen.
                                finalPopSizes[i] = 0.0; }
                            else {
                                finalPopSizes[i] = finalPopSizes[i] - oldGT;
                            }}}
                }
            }
        }
        //Now combine the new mutants with the deterministic growth
        ElementAddInPlace(finalPopSizes, newMuts);
        double finalPopSize = std::accumulate(finalPopSizes.begin(), finalPopSizes.end(), 0.0); // This should be very similar to what we calculated before mutations came into play
        //Now sample in to the next generation after the transfer.
        double expectationFactor = N0 / finalPopSize;
        for (int i = 0; i < finalPopSizes.size(); i++) {
            auto cPopSize = finalPopSizes[i];
            if (cPopSize > 0)
            {
                double expectation = cPopSize * expectationFactor;
                PopSizes[i] = RPOIS(cPopSize);
            } else{
              PopSizes[i] = 0.0;
            }
        }
}

}
