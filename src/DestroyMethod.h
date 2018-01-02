//
// Created by alberto on 27/06/17.
//

#ifndef ML_PALNS_DESTROYMETHOD_H
#define ML_PALNS_DESTROYMETHOD_H

#include <random>
#include <memory>

namespace mlpalns {
    // Destroy method
    template<class Solution>
    struct DestroyMethod {
        DestroyMethod() = default;

        virtual void destroy_solution(Solution& sol, std::mt19937& mt) = 0;

        virtual std::unique_ptr<DestroyMethod<Solution>> clone() const = 0;

        virtual ~DestroyMethod() = default;
    };
}

#endif // ML_PALNS_DESTROYMETHOD_H
