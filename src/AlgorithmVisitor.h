//
// Created by alberto on 28/06/17.
//

#ifndef ML_PALNS_ALGORITHMVISITOR_H
#define ML_PALNS_ALGORITHMVISITOR_H

#include "AlgorithmStatus.h"
#include "DestroyMethod.h"
#include "RepairMethod.h"
#include <vector>
#include <string>

namespace mlpalns {
    /**
     * This class implements a Visitor for the algorithm, to be called at
     * specific points during the solution process.
     *
     * @tparam Solution The solution type.
     */
    template<class Solution>
    struct AlgorithmVisitor {
        /**
         * This method is called at the start of the algorithm, before the first iteration.
         * The vectors passed to this method can be indexed via AlgorithmStatus::destroy_method_id
         * and AlgorithmStatus::repair_method_id.
         *
         * @param destroy A vector with the destroy methods.
         * @param repair  A vector with the repair methods.
         * @param destroy_methods_desc A vector with string descriptions of the destroy methods.
         * @param repair_methods_desc  A vector with string descriptions of the repair methods.
         */
        virtual void on_algorithm_start(std::vector<DestroyMethod<Solution>*>& destroy,
                                        std::vector<RepairMethod<Solution>*>& repair,
                                        const std::vector<std::string>& destroy_methods_desc,
                                        const std::vector<std::string>& repair_methods_desc) = 0;

        /**
         * This method is called when the prerun (calibration run) ends.
         *
         * @param destroy A vector with the destroy methods.
         * @param repair  A vector with the repair methods.
         */
        virtual void on_prerun_end(std::vector<DestroyMethod<Solution>*>& destroy,
                                   std::vector<RepairMethod<Solution>*>& repair) = 0;

        /**
         * This method is called at the end of every iteration.
         *
         * @param alg_status    An AlgorithmStatus object containing info on the current state of the algorithm.
         */
        virtual void on_iteration_end(AlgorithmStatus<Solution>& alg_status) = 0;

        /**
         * This method is called when many iterations have passed, without any improvement
         * on the objective value of the best known solution. In this case, it is assumed that
         * the user might want to change in the repair and destroy methods.
         * How many iterations need to pass is determined by Params::iters_without_improvement_alarm.
         *
         * @param destroy A vector with the destroy methods.
         * @param repair  A vector with the repair methods.
         */
        virtual void on_many_iters_without_improvement(std::vector<DestroyMethod<Solution>*>& destroy,
                                                       std::vector<RepairMethod<Solution>*>& repair) = 0;

        /**
         * Virtual destructor.
         */
        virtual ~AlgorithmVisitor() = default;
    };
}

#endif // ML_PALNS_ALGORITHMVISITOR_H
