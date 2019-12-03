//
// Created by alberto on 28/06/17.
//

#ifndef ML_PALNS_ALGORITHMSTATUS_H
#define ML_PALNS_ALGORITHMSTATUS_H

#include <cstdint>
#include <memory>
#include <vector>

#include "Parameters.h"

namespace mlpalns {
    template<class Solution>
    struct AlgorithmStatus {
        // Params object
        const Parameters& params;

        // New solution produced at this iter
        Solution& new_solution;

        // Global best solution
        Solution& best_solution;

        // Iteration number
        std::uint32_t iter_number;

        // Id of the destroy method chosen at this iter
        std::size_t destroy_method_id;

        // Id of the repair method chosen at this iter
        std::size_t repair_method_id;

        // Whether the new sol was accepted at this iter
        bool accepted;

        // Whether the new sol was the new global best sol
        bool new_best;

        AlgorithmStatus(const Parameters& params, Solution& new_solution, Solution& best_solution)
            : params{params}, new_solution{new_solution}, best_solution{best_solution} {}
    };
} // namespace mlpalns

#endif // ML_PALNS_ALGORITHMSTATUS_H
