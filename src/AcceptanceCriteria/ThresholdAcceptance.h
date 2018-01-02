#ifndef THRESHOLD_ACCEPTING_H
#define THRESHOLD_ACCEPTING_H

#include "../AcceptanceCriterion.h"

namespace mlpalns {
    /*! \brief This class models the Threshold Accepting acceptance criterion for ALNS */
    template<typename Solution>
    struct ThresholdAcceptance : public AcceptanceCriterion<Solution> {
        /*! The current gap (threshold) */
        double current_threshold;

        /*! The amount by which the threshold is decreased, if the decrease is linear */
        double threshold_decrease;

        /*! The coefficient by which the threshold is multiplied, if the decrease is exponential */
        double threshold_exp_coeff;

        /*! Basic constructor */
        ThresholdAcceptance(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

        /*! This method initialise method-specific parameters */
        void initialise() override;

        /*! This method is used to calibrate the parameters, given a starting solution's data */
        void calibrate(const CalibrationData& d) override;

        /*! This method provides some info about the current status of the acceptance criterion, to be printed out during the computation. In particular, we
         * will print out the temperature! */
        std::string get_print_str() const override;

        /*! This method updates the acceptance criterion's parameters, based on running info
         *
         *  @param iter_number is the current iteration number
         *  @param best_obj is the value of the current best solution
         */
        void update_parameters(std::uint32_t iter_number, double best_obj) override;

        /*! This method returns true iff the solution should be accepted according to the acceptance criterion
         *
         *  @param best_obj is the objective value of the best solution (globally) found so far
         *  @param current_obj is the objective value of the current solution
         *  @param new_obj is the objective value of the new solution (which we must determine wether to accept or not)
         *  @param eps is an epsilon to be used in numerical comaprisons
         *  @param new_sol is a reference to the new solution found
         *  @param mt is a reference to a random number generator
         */
        bool should_accept(double best_obj, double current_obj, double new_obj, double eps, Solution& new_sol, std::mt19937& mt) override;
    };

    template<typename Solution>
    void ThresholdAcceptance<Solution>::initialise() {
        current_threshold = 1.0;
        threshold_decrease = 0.0;
        threshold_exp_coeff = 1.0;
    }

    template<typename Solution>
    void ThresholdAcceptance<Solution>::calibrate(const CalibrationData&) {
        current_threshold = this->params.ta_params.start_threshold;
        threshold_decrease =
            (this->params.ta_params.start_threshold - this->params.ta_params.end_threshold) / (this->params.max_iters - this->params.prerun_iters);
        threshold_exp_coeff =
            std::pow(this->params.ta_params.start_threshold / this->params.ta_params.end_threshold, 1.0 / (this->params.max_iters - this->params.prerun_iters));
    }

    template<typename Solution>
    std::string ThresholdAcceptance<Solution>::get_print_str() const {
        return ("threshold: " + std::to_string(static_cast<long double>(current_threshold)));
    }

    template<typename Solution>
    void ThresholdAcceptance<Solution>::update_parameters(std::uint32_t, double) {
        if(this->params.ta_params.threshold_decrease_is_linear) {
            current_threshold -= threshold_decrease;
        } else {
            current_threshold *= threshold_exp_coeff;
        }
    }

    template<typename Solution>
    bool ThresholdAcceptance<Solution>::should_accept(double, double current_obj, double new_obj, double eps, Solution&, std::mt19937&) {
        double gap = (new_obj - current_obj) / new_obj;
        return (gap < current_threshold - eps);
    }
}

#endif