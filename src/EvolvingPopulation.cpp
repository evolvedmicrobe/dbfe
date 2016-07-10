#include "Classes.h"
#include <numeric>
#include <math>

namespace dbfe {


   inline dvector Element_A_Times_E_ToThe_RT(const dvector& A, const dvector&  r, const double t)
   {
        if (A.size() != r.size())
        { Rpp::stop("Cannot Multiply vectors of Unequal Length"); }
        dvector toR(A.length()); 
        for (int i = 0; i < A.size(); i++)
        {
            toR[i] = A[i] * exp(r[i] * t) ;
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
        if (x.length() != y.length())
        { Rcpp::stop("Cannot Subtract and Sum Vectors of Unequal Length"); }
        for (int i = 0; i < x.length(); i++)
        {
            x[i] += y[i];
        }
    }

	EvolvingPopulation::EvolvingPopulation(spDiscretizedDFE dfe, const PopulationSize& popSize) :
	 dfe(dfe),
	 PopSizes(dfe->NumberOfClassesIncludingNeutral()),
	 popRelativeFitnesses(dfe->MidPoints),
	 MutCounter(dfe)
	  {	
          LogrelativePopIncrease = std::log(popSize.NF / popSize.N0);
          InvRelativePopIncrease = popSize.N0 / popSize.NF;
          PopSizes[0] = popSize.N0;
          N0 = popSize.N0;
          NF = popSize.NF;

	  }

	EvolvingPopulation::GrowOneCycle() {

        double curPopSize = std::accumulate(PopSizes.begin(), PopSizes.end(), 0.0);
        dvector freqs(curPopSize.size());
        std::transform(PopSizes.begin(), PopSizes.end(), freqs.begin(), [curPopSize](const double a) {return a / curPopSize;} );
        // Get average fitness
        double meanFitnessStart = std::inner_product(freqs.begin(), freqs.end(), popRelativeFitnesses.begin(), 0.0);

        //how much growth occurs, approximate by how long the mean fitness would take to grow
        double growthTime = log(NF / curPopSize) / meanFitnessStart;
        //First Handle Mutations
        dvector newMuts(PopSizes.size());
        //dvector finalPopSize =  std::transform(PopSizes.begin(), PopSizes.end(), popRelativeFitnesses.begin(), [growthTime] )            
        dvector finalPopSize = Element_A_Times_E_ToThe_RT(PopSizes, popRelativeFitnesses, growthTime);
        dvector NFminusN0akaTimeForPoisson  = ElementSubtract(FinalPopSize, PopSizes);//also the time for poisson.
        //Sum of N_f- N_0, note approximate as if a mutation occurs that one guy that mutates will no longer contribute to this time
        TotalTime += std::accumulate(NFminusN0akaTimeForPoisson.begin(), NFminusN0akaTimeForPoisson.end(), 0.0);
        
        dvector& mutRates = dfe.MutationRates;
        for (int i = 0; i < PopSizes.size(); i++)
        {
            //Only work with populations that exist to mutate
            double popSize = PopSizes[i];
            if (popSize > 0)
            {                 
                //New algorithm, for each class we will go through and mutate
                //the expected number of each fitness class
                //TODO: Decide later if in the case of very few classes (mean muts =0) we should not do one poisson for all classes and 
                //assign after the fact. Note small poissons are ~3X more expensive than single unifrom with the MersenneTwister RNG
                for(int j = 0; j < mutRates.Length; j++) {
                    //Determine the mean number of mutations
                    double meanMutForClass = mutRates[j] * popSize;                    
                    //Sample poisson random
                    int mutNumber = RPOIS(meanMutForClass);
                    int w = j + 1;                        
                    //Get the maximum time on the rescaled valued, and subtract one as it ranges from 1 to this high value, so this
                    //will be the multiplication factor
                    MutCounter.CountOfEachMutation[w] += mutNumber;
                    //If higher, add to population
                    if (w > i)//Not the *i* index includes the zero class, while the *j* index does not, so corresponds to i+1
                    {
                        //TODO: Why is 1 here?  Did I do that to account for the need of one division to happen?
                        //time seems to be on scale [1,GrowthTime]
                        //double maxRescaled = Math.Exp(GrowthTime * dfe.MidPoints[i]) - 1.0;                        
                        double maxRescaled = exp(GrowthTime * dfe.MidPoints[i]);                        
                        
                        //Now generate a mutation for each class
                        for (int k = 0; k < mutNumber; k++)
                        {
                            //double uniformTime = 1.0 + RandomVariateGenerator.NextDouble() * maxRescaled;
                            double uniformTime = RUNIF * maxRescaled;
                            //Convert back
                            double actualTime = log(uniformTime) / dfe.MidPoints[i];
                            //Find out how much the mutant grew, then add it to the new group
                            //and subtract from the old population group the growth it would have contributed there
                            //had it not mutated
                            double gr = dfe.MidPoints[w];
                            double GT = exp(gr * (GrowthTime - actualTime));
                            newMuts[w] += GT;
                            double old_gr = dfe.MidPoints[i];
                            double oldGT = exp(old_gr * (GrowthTime - actualTime));
                            if (oldGT > FinalPopSize[i]) //if the growth is higher due to the approximation, subtract it all, should rarely happen.
                            {
                                FinalPopSize[i] = 0.0; }
                            else {
                                FinalPopSize[i] = FinalPopSize[i] - oldGT;
                            }}}}}}
            //Now combine the new mutants with the deterministic growth
            ElementAddInPlace(newMuts, FinalPopSize); //do as two steps to avoid array allocation
            // TODO: Verify this copy works okay, perhaps std::move would be better?
            PopSizes = newMuts;
            //Now sample in to the next generation after the transfer.
            double expectationFactor = N0 / std::accumulate(PopSizes.begin(), Sum();
            for (int i = 0; i < PopSizes.Length; i++) {
                var cPopSize=PopSizes[i];
                if (cPopSize> 0)
                {
                    double expectation = cPopSize * expectationFactor;
                    PopSizes[i] = RandomVariateGenerator.PoissonSample(cPopSize);
                }}
	}

}