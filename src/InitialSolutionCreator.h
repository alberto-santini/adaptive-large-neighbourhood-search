//
// Created by alberto on 27/06/17.
//

#ifndef ML_PALNS_INITIALSOLUTIONCREATOR_H
#define ML_PALNS_INITIALSOLUTIONCREATOR_H

#include <random>

namespace mlpalns {
    // Constructive heuristic to get an initial solution
    template<class Solution, class ProblemInstance>
    struct InitialSolutionCreator {
        virtual Solution create_initial_solution(const ProblemInstance& instance, std::mt19937& mt) = 0;
    };
}

#endif // ML_PALNS_INITIALSOLUTIONCREATOR_H
