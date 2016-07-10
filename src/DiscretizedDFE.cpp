#include "DiscretizedDFE.h"


using namespace std;

namespace dbfe {

    DiscretizedDFE::DiscretizedDFE (double min, double max, int numberBins) :
    min(min),
    step(0),
    MidPoints(numberBins + 1),
    _pClassProbabilities(numberBins),
    _pMutationRates(numberBins),
    cumProbs(numberBins)
    {
    	this->max = (1 + max) * ANCESTRALGROWTHRATE;
 		double interval = (max - min) / numberBins;
 		MidPoints[1] = min + interval / 2.0;
 		//Make a list of relative fitness improvements
	    for (int i = 2; i < numberBins + 1; i++) {
			MidPoints [i] = MidPoints [i - 1] + interval;
		}
		// Now just get the actual growth rates
		for (int i = 0; i < numberBins + 1; i++) {
			MidPoints [i] = ANCESTRALGROWTHRATE + MidPoints [i] * ANCESTRALGROWTHRATE;
		}       
    }

	DiscretizedDFE::DiscretizedDFE (double max, int numberBins) :
	MidPoints(numberBins + 1),
	ClassProbabilities(numberBins),
	_pMutationRates(numberBins),
	cumProbs(numberBins),
	min(0.0)
		{
			this.max = (1 + max) * ANCESTRALGROWTHRATE;
			double interval = max / NumberBins;
			MidPoints [1] = interval / 2.0;
			//Make a list of relative fitness improvements
			for (int i = 2; i < NumberBins + 1; i++) {
				MidPoints [i] = MidPoints [i - 1] + interval;
			}
			//Now just get the actual growth rates
			for (int i = 0; i < NumberBins + 1; i++) {
				MidPoints [i] = ANCESTRALGROWTHRATE + MidPoints [i] * ANCESTRALGROWTHRATE;
			}
		}

	void DiscretizedDFE::CreateCumulativeProbs () {
		cumProbs[0] = _pClassProbabilities[0];
        for(int i=1; i<_pClassProbabilities.Length; i++)
        {
            cumProbs[i] = _pClassProbabilities[i] + cumProbs[i - 1];
        }
	}

	void DiscretizedDFE::SetProb(const dvector& probs, double mutationRate);
        {
            if (probs.size() != _pClassProbabilities.size())
            {
                Rcpp::stop("Probability arrays were of unequal size");
            }
            for (int i = 0; i < probs.size(); i++)
            {
                _pClassProbabilities[i] = probs[i];
                _pMutationRates[i] = probs[i] * mutationRate;
            }
            CreateCumulativeProbs();
        }

    int DiscretizedDFE::AssignFitnessToBin (double W)
        {
            if (W < min) {
                Rcpp::stop("Can't assign below the min!");
            }
            W = W * ANCESTRALGROWTHRATE;
            double dif = abs(W - MidPoints[i])
            for(int i = 1; i < MidPoints.size(); i++) {
            	double new_dif = abs(W - MidPoints[i]);
            	if (new_dif > dif) {
            		return i - 1;
            	}
            }
            return MidPoints.size() - 1;
        }
    void DiscretizedDFE::UpdateWithNewSamples (const std::vector<ObservedWell>& AugmentedData)
        {
            ///Sample the rate for each class, then update the cumulative probabilities
            double rateTotal = 0.0;
            double totalTime=  std::accumulate(AugmentedData.begin(), AugmentedData.end(), 0.0,
            	[](const ObservedWell& ow, double b) {
            		return b + ow.AmountOfTimeLastRun;});
            //Sample a new mu
            //Gamma is parameterized with shape and rate (rate = 1/scale)
            for (int i = 0; i < _pMutationRates.size(); i++)
            {
            	// AugmentedData.Select(x => x.MutCounter.CountOfEachMutation[i - 1]).Sum();                
                int totalMutations = std::accumulate(AugmentedData.begin(), AugmentedData.end(), 0, 
                	[i](const ObservedWell& ow, int b) {
                		return b + ow.MutCounter.CountOfEachMutation[i - 1];
                	});
                // Should take a shape and scale (which is the inverse of the rate)
                double newRate= R::rgamma((priorAlpha + static_cast<double>(totalMutations)), 1.0 / (priorBeta + totalTime));
                rateTotal += newRate;
                _pMutationRates[i] = newRate;  
            }
            // ClassProbabilities = _pMutationRates.ElementDivide(rateTotal);
            std::transform(_pMutationRates.begin(), _pMutationRates.end(), _pMutationRates.begin(),
            	[rateTotal](const double x) {return x / rateTotal;});
            
            CreateCumulativeProbs ();
            PointMass = false;
        }

// Going to shift this to R.
     void DiscretizedDFE::InitializeWithObservedData(const std::vector<ObservedWell>& Data)
        {
            //A mutation escapes loss at probability ~2s, and it rises to high frequency with probability 
            //Ne = N0*growthrate*t or N0*ln(2)*(Gen. Between Transfers).
            //time to hit X% frequency is: T= ln((X/10)*(Ne-1))/s

            //so here is the strategy, ignoring clonal interferance, the expected number of each class present is equal to the 
            //number that escape drift each generation times the frequency they achieve after X amount of time, so will set these two equal

            //number that appear each generation = Ne*mu*2s
            //end frequency = 1/ (1+exp(-s*t)*(Ne-1)
            //so expectated frequency = Sum_(over generations) Ne*mu*2s*(1/(1+exp(-s*t)*(Ne-1))
            //and this can be solved for the initial rate.


//            var totTransfers = Data.Select(x => x.TotalTransfers).Distinct();
//            if (totTransfers.Count() > 1)
//            {
//                throw new Exception("Code needs to be changed to handle different number of transfers.");
//            }
//            //otherwise just go by variables
//            var popSizes = Data.GroupBy(y => y.PopSize).Select(x => new { PopSize = x.Key, Obs = x.ToList() }).ToList();
//            for(int i=1;i<NumberOfClassesIncludingNeutral;i++)
//            {
//                double s = (MidPoints[i] / MidPoints[0]) - 1;
//                //estimate rate for each pop size
//                double muEst = 0;// new double[popSizes.Count];
//                List<double> estimates = new List<double>();
//                foreach(var curPop in popSizes)
//                {
//                    var ps=curPop.PopSize;
//                    var Ne = ps.Ne;
//                    var totalGens = curPop.Obs[0].TotalGenerations;
//                    var ne2s = Ne * 2 * s;
//                    double[] expectations = new double[(int)totalGens];
//                    var expectedFreqSansMu=ne2s*Enumerable.Range(0,(int)totalGens).Select(t=>  1 / (1+Math.Exp(-s*(totalGens-t))*(Ne-1))).Sum();
//                    var obsFreq=curPop.Obs.Where(x=>x.BinClass==i).Count()/(double)curPop.Obs.Count;
//                    var curEst = obsFreq / expectedFreqSansMu;
//                    estimates.Add(curEst);
//                    muEst += curEst / popSizes.Count;//get average by summing x*(1/N)
//                }
//                estimates.Add(1e-20);
//                _pMutationRates[i - 1] = estimates.Max();
//            }
        }
        void DiscretizedDFE::SetDFEasPointMass (int IndexToSet)
        {
            std::fill(_pClassProbabilities.begin(), _pClassProbabilities.end(), 0.0);
            _pClassProbabilities [IndexToSet] = 1.0;
            CreateCumulativeProbs ();
            PointMass = true;
            PointMassGroup = IndexToSet + 1;
        }       

}