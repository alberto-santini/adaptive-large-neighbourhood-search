//
// Created by alberto on 27/06/17.
//

#ifndef ML_PALNS_REPAIRMETHOD_H
#define ML_PALNS_REPAIRMETHOD_H

#include <random>
#include <memory>

namespace mlpalns {
    // Repair method
    template<class Solution>
    struct RepairMethod {
        RepairMethod() = default;

        virtual void repair_solution(Solution& sol, std::mt19937& mt) = 0;

        virtual std::unique_ptr<RepairMethod<Solution>> clone() const = 0;

        virtual ~RepairMethod() = default;
    };
}

#endif // ML_PALNS_REPAIRMETHOD_H
