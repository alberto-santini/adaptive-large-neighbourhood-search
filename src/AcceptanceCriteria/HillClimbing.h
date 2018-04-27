#ifndef HILL_CLIMBING_H
#define HILL_CLIMBING_H

#include "../AcceptanceCriterion.h"

namespace mlpalns {
    /*! \brief This class models the Hill Climbing acceptance criterion for ALNS */
    template<typename Solution>
    struct HillClimbing : public AcceptanceCriterion<Solution> {
        /*! Basic constructor */
        HillClimbing(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void HillClimbing<Solution>::initialise() {}

    template<typename Solution>
    void HillClimbing<Solution>::calibrate(const CalibrationData&) {}

    template<typename Solution>
    std::string HillClimbing<Solution>::get_print_str() const {
        return std::string();
    }

    template<typename Solution>
    void HillClimbing<Solution>::update_parameters(std::uint32_t, double, double) {}

    template<typename Solution>
    bool HillClimbing<Solution>::should_accept(double, double current_obj, double new_obj, double eps, Solution&, std::mt19937&) {
        return (new_obj < current_obj - eps);
    }
}

#endif