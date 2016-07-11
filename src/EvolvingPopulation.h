#ifndef EvolvingPopulation_h
#define EvolvingPopulation_h

#include "Common.h"
#include <memory>

#include "DiscretizedDFE.h"
#include "MutationCounter.h"
#include "ObservedWell.h"

namespace dbfe {

    typedef std::shared_ptr<DiscretizedDFE> spDiscretizedDFE;

    class EvolvingPopulation   {
    public:
        std::vector<double> PopSizes;
        std::vector<double> popRelativeFitnesses;
        double N0;
        double NF;
        double LogrelativePopIncrease;
        double InvRelativePopIncrease;
        double TotalTime;
        spDiscretizedDFE dfe;
        MutationCounter MutCounter;

        EvolvingPopulation(spDiscretizedDFE dfe, const PopulationSize& popSize);
        EvolvingPopulation(spDiscretizedDFE dfe, ObservedWell well) :
        EvolvingPopulation(dfe, *well.PopSize) {};

        void GrowOneCycle();
        void GrowOneCycleLabeled();
        int SamplePopulation();
    };
}
#endif /* EvolvingPopulation_h */

