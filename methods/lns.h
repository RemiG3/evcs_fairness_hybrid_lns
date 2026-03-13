#ifndef LNS_H
#define LNS_H

#include "utils.h"


class LargeNeighboroodSearch {
public:
    double initialTemperature;
    double finalTemperature;
    int maxTrials;
    int maxGenerated;
    int lastBestObjective;

    bool save_objective_evolution;
    vector<double> objectiveEvolution;
    int seed;

    LargeNeighboroodSearch(double initialTemp = 10.0, double finalTemp = 0.00001, int trials = 1000, int generated = 100, bool save_obj_evol = false, int seed = -1);
    void initializeSearch(AnnealingState &state, ProblemParams const &params);
    void constructiveOperator(AnnealingState& state, const ProblemParams& params);
    void destructiveOperator(AnnealingState& state, const ProblemParams& params);
    void searchFunction(AnnealingState& state, const ProblemParams& params);
    double objective(AnnealingState& state, const ProblemParams& params);
    void optimize(ProblemParams const &params, AnnealingState &state, double timeout);
    void saveObjectiveEvolutionToFile(const std::string &filepath);
};


void repairingSolution(AnnealingState& state, const ProblemParams& params);



#endif // LNS_H 