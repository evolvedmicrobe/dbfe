//
//  DiscretizedDFE.h
//  dbfe
//
//  Created by Nigel Delaney on 7/10/16.
//
//

#ifndef DiscretizedDFE_h
#define DiscretizedDFE_h

#include "Common.h"

namespace dbfe
{

const double ANCESTRALGROWTHRATE = 0.693147180559945;

class DiscretizedDFE {

    protected:
        double max, min, step;
        std::vector<double> _pClassProbabilities;
        /// <summary>
        /// The rates in each class
        /// </summary>
        std::vector<double> _pMutationRates;

        /// <summary>
        /// 1 to N
        /// </summary>
        dvector cumProbs;

        void CreateCumulativeProbs ();

    private:
        bool PointMass = false;
        int PointMassGroup = 0;

    public:
        /// <summary>
        /// 0 to N, GROWTH RATES! Not fitness per say
        /// </summary>
        std::vector<double> MidPoints;
        double priorAlpha = 1.0;
        double priorBeta = 1.0;
        /// <summary>
        /// The rate obtained by summing the individual rates.
        /// </summary>
        double BeneficialMutationRate() { return std::accumulate(_pMutationRates.begin(), _pMutationRates.end(), 0.0);};

        /// <summary>
        /// The rate at which mutations appear per unit time per individual.
        /// This array does NOT have the neutral class.
        /// </summary>
        std::vector<double>& MutationRates() { return _pMutationRates;  };

        /// <summary>
        /// 1 to N
        /// </summary>
        std::vector<double>& ClassProbabilities() {return _pClassProbabilities;} // getter
        void ClassProbabilities(std::vector<double> values); //setter

        size_t NumberOfClassesIncludingNeutral() const {return MidPoints.size();}

        double GetGrowthRateForBin (int binNumber) { return MidPoints [binNumber + 1];}


        static double ConvertGrowthRateToSelectiveCoefficient(double gr)
        {
            return (gr / ANCESTRALGROWTHRATE) - 1;
        }

        dvector SelectiveCoefficients() {
            dvector coefs(MidPoints.size());
            std::transform(MidPoints.begin(), MidPoints.end(), coefs.begin(), ConvertGrowthRateToSelectiveCoefficient);
            return coefs;
        };


        /// <summary>
        /// Returns a random bin for a mutation based on the
        /// current parameter values, BIN NUMBER IS in 0 to 1 scheme
        /// </summary>
        /// <returns></returns>
        int GetRandomBinAssignment () const
        {

            double d = RUNIF;
            for (int i = 0; i < cumProbs.size(); i++) {
                if (d < cumProbs [i])
                    return i + 1;
            }
            Rcpp::stop("Problem with getting random bin");
            return 0; // Never happens, but gets rid of warning.
        };

        /// <summary>
        /// Overwrite the values in the existing array instead of copying and allocating.
        /// </summary>
        /// <param name="probs"></param>
        void SetProb(const dvector& probs, const double mutationRate);

        //Fitness is s, where s is (r0+s)/r0;
        /// <summary>
        /// Assign a fitness value
        /// </summary>
        /// <param name="W">1+s, to be multiplied by ancestral growth rate</param>
        /// <returns></returns>
        int AssignFitnessToBin (double W) const;


        void UpdateWithNewSamples (const std::vector<ObservedWell>& AugmentedData);

        /// <summary>
        /// Tries to guess good initial parameters based on some observed values.
        /// </summary>
        /// <param name="Data"></param>
        void InitializeWithObservedData(const std::vector<ObservedWell>& Data);

        void SetDFEasPointMass (int IndexToSet);


        // CONSTRUCTORS
        DiscretizedDFE (double min, double max, int NumberBins);

        DiscretizedDFE (double max, int NumberBins);

    };
}; // end namespace

#endif /* DiscretizedDFE_h */
