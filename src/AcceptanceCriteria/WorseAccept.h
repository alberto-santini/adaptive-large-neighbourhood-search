#ifndef WORSE_ACCEPT_H
#define WORSE_ACCEPT_H

#include "../AcceptanceCriterion.h"

namespace mlpalns {
    /*! \brief This class models the Worse Accept acceptance criterion for ALNS */
    template<typename Solution>
    struct WorseAccept : public AcceptanceCriterion<Solution> {
        /*! The current chance of accepting a worsening solution */
        double current_prob;

        /*! Are we decreasing the probability linearly or exponentially? */
        bool prob_decrease_is_linear;

        /*! If the probability decreases linearly, that's the amount by which it decreases at each iteration */
        double prob_decrease;

        /*! Basic constructor */
        WorseAccept(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void WorseAccept<Solution>::initialise() {
        current_prob = 0.1;
        prob_decrease = 0.0;
        prob_decrease_is_linear = true;
    }

    template<typename Solution>
    void WorseAccept<Solution>::calibrate(const CalibrationData&) {
        current_prob = this->params.wa_params.start_prob;
        prob_decrease_is_linear = this->params.wa_params.prob_decrease_is_linear;
        prob_decrease = (current_prob - this->params.wa_params.end_prob) / (this->params.max_iters - this->params.prerun_iters);
    }

    template<typename Solution>
    std::string WorseAccept<Solution>::get_print_str() const {
        return ("worsening probability: " + std::to_string(static_cast<long double>(current_prob)));
    }

    template<typename Solution>
    void WorseAccept<Solution>::update_parameters(std::uint32_t iter_number, double) {
        if(prob_decrease_is_linear) {
            current_prob -= prob_decrease;
        } else {
            double S = this->params.wa_params.start_prob;
            double E = this->params.wa_params.end_prob;
            auto N = this->params.max_iters - this->params.prerun_iters;
            double lambda = log(S / E) / N;
            auto n = iter_number - this->params.prerun_iters;
            current_prob = S * exp(-lambda * n);
        }
    }

    template<typename Solution>
    bool WorseAccept<Solution>::should_accept(double, double current_obj, double new_obj, double eps, Solution&, std::mt19937& mt) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return (new_obj < current_obj - eps || dist(mt) < current_prob - eps);
    }
}

#endif