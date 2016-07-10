



    typedef std::shared_ptr<DiscretizedDFE> spDiscreitzedDFE;

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
        EvolvingPopulation(spDiscretizedDFE dfe, ObservedWell well) :this(dfe,well.PopSize);
        
        void GrowOneCycle();
        void GrowOneCycleLabeled();
        int SamplePopulation();       
    };


    class MutationCounter {
        public: 
            int CountOfEachMutation[];
            MutationCounter(spDiscretizedDFE dfe);

    }

}
