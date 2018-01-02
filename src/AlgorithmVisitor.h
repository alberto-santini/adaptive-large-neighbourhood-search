//
// Created by alberto on 28/06/17.
//

#ifndef ML_PALNS_ALGORITHMVISITOR_H
#define ML_PALNS_ALGORITHMVISITOR_H

#include "AlgorithmStatus.h"

namespace mlpalns {
    template<class Solution>
    struct AlgorithmVisitor {
        AlgorithmVisitor() = default;
        
        virtual void on_iteration_end(AlgorithmStatus<Solution>& alg_status) = 0;
        
        virtual ~AlgorithmVisitor() = default;
    };
}

#endif // ML_PALNS_ALGORITHMVISITOR_H
