#ifndef LATE_ACCEPTANCE_HILL_CLIMBING_H
#define LATE_ACCEPTANCE_HILL_CLIMBING_H

#include "../AcceptanceCriterion.h"
#include <limits>
#include <sstream>

namespace mlpalns {
    /*! \brief This class models the Late Acceptance Hill Climbing acceptance criterion for ALNS */
    template<typename Solution>
    struct LateAcceptanceHillClimbing : public AcceptanceCriterion<Solution> {
        /*! The list of saved solutions */
        std::vector<double> sol_list;

        /*! Current selected element in the list */
        std::vector<double>::size_type sol_index;

        /*! Basic constructor */
        LateAcceptanceHillClimbing(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void LateAcceptanceHillClimbing<Solution>::initialise() {
        sol_list = std::vector<double>(this->params.lahc_params.list_size, std::numeric_limits<double>::max());
        sol_index = 0;
    }

    template<typename Solution>
    void LateAcceptanceHillClimbing<Solution>::calibrate(const CalibrationData& d) {
        sol_list = std::vector<double>(this->params.lahc_params.list_size, d.objective_value);
        sol_index = 0;
    }

    template<typename Solution>
    std::string LateAcceptanceHillClimbing<Solution>::get_print_str() const {
        double current_late_acceptance_value = sol_list[sol_index];
        std::stringstream val;

        val << "compare to value[" << sol_index << "] = ";

        if(current_late_acceptance_value == std::numeric_limits<double>::max()) {
            val << "infinite";
        } else {
            val << current_late_acceptance_value;
        }

        return val.str();
    }

    template<typename Solution>
    void LateAcceptanceHillClimbing<Solution>::update_parameters(std::uint32_t, double, double) {
        sol_index++;

        if(sol_index == sol_list.size()) {
            sol_index = 0;
        }
    }

    template<typename Solution>
    bool LateAcceptanceHillClimbing<Solution>::should_accept(double, double, double new_obj, double eps, Solution&, std::mt19937&) {
        bool accepted = false;
        long int current_index = static_cast<long int>(sol_index) - 1;

        if(current_index < 0) {
            current_index = sol_list.size() - 1;
        }

        if((new_obj < sol_list[sol_index] - eps) || (this->params.lahc_params.accept_non_worsening && new_obj < sol_list[current_index] - eps)) {
            accepted = true;
        }

        if(accepted) {
            sol_list[sol_index] = new_obj;
        }

        return accepted;
    }
} // namespace mlpalns

#endif