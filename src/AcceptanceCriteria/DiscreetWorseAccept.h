#ifndef DISCREET_WORSE_ACCEPT_H
#define DISCREET_WORSE_ACCEPT_H

#include "../AcceptanceCriterion.h"

namespace mlpalns {
    /*! \brief This class models the Conservative Worse Accept acceptance criterion for ALNS */
    template<typename Solution>
    struct DiscreetWorseAccept : public AcceptanceCriterion<Solution> {
        /*! The current chance of accepting a worsening solution */
        double current_prob;

        /*! Number of consecutive rejections of new vs current */
        std::uint32_t consecutive_reject;

        const double PROBABILITY_ZERO;

        /*! Basic constructor */
        DiscreetWorseAccept(Parameters& par) : AcceptanceCriterion<Solution>(par), PROBABILITY_ZERO(10e-7) {}

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
         *  @param elapsed_time is the current elapsed time
         *  @param best_obj is the value of the current best solution
         */
        void update_parameters(std::uint32_t iter_number, double elapsed_time, double best_obj) override;

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
    void DiscreetWorseAccept<Solution>::initialise() {
        current_prob = 0.1;
        consecutive_reject = 0;
    }

    template<typename Solution>
    void DiscreetWorseAccept<Solution>::calibrate(const CalibrationData&) {
        current_prob = 0.0;
        consecutive_reject = 0;
    }

    template<typename Solution>
    std::string DiscreetWorseAccept<Solution>::get_print_str() const {
        return ("worsening probability: " + std::to_string(static_cast<long double>(current_prob)));
    }

    template<typename Solution>
    void DiscreetWorseAccept<Solution>::update_parameters(std::uint32_t, double, double) {
        auto N = this->params.dwa_params.consecutive_rejects_for_100p;

        if(this->params.dwa_params.prob_increase_is_linear) {
            current_prob = (double)consecutive_reject / (double)N;
        } else {
            double lambda = log(1 / PROBABILITY_ZERO) / (double)N;
            current_prob = PROBABILITY_ZERO * exp(lambda * consecutive_reject);
        }
    }

    template<typename Solution>
    bool DiscreetWorseAccept<Solution>::should_accept(double, double current_obj, double new_obj, double eps, Solution&, std::mt19937& mt) {
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        bool accept = (new_obj < current_obj - eps || dist(mt) < current_prob - eps);

        if(!accept) {
            ++consecutive_reject;

            if(consecutive_reject > this->params.dwa_params.consecutive_rejects_for_100p) {
                consecutive_reject = this->params.dwa_params.consecutive_rejects_for_100p;
            }
        } else {
            consecutive_reject = 0;
        }

        return accept;
    }
} // namespace mlpalns

#endif