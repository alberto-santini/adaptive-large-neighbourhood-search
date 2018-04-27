#ifndef NL_GREAT_DELUGE_H
#define NL_GREAT_DELUGE_H

#include "../AcceptanceCriterion.h"
#include <cmath>
#include <limits>

namespace mlpalns {
    /*! \brief This class models the Great Deluge acceptance criterion for ALNS */
    template<typename Solution>
    struct NLGreatDeluge : public AcceptanceCriterion<Solution> {
        /*! The water level. We accept only "underwater" solutions! */
        double water_level;

        /*! Objective value of the current solution */
        double current_obj;

        /*! Objective value of the best solution */
        double best_obj;

        /*! Basic constructor */
        NLGreatDeluge(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void NLGreatDeluge<Solution>::initialise() {
        water_level = std::numeric_limits<double>::max();
        current_obj = 0;
        best_obj = 0;
    }

    template<typename Solution>
    void NLGreatDeluge<Solution>::calibrate(const CalibrationData& d) {
        water_level = d.objective_value * this->params.nlgd_params.nl_initial_water_level_ratio;
        current_obj = d.objective_value;
        best_obj = d.objective_value;
    }

    template<typename Solution>
    std::string NLGreatDeluge<Solution>::get_print_str() const {
        std::string wl = "";

        if(water_level == std::numeric_limits<double>::max()) {
            wl = "infinite";
        } else {
            wl = std::to_string(static_cast<long double>(water_level));
        }
        return ("water level: " + wl);
    }

    template<typename Solution>
    void NLGreatDeluge<Solution>::update_parameters(std::uint32_t, double, double best_obj) {
        double gap = (water_level - current_obj) / water_level;

        if(gap < this->params.nlgd_params.nl_gap_to_increase_water_level) {
            water_level += std::abs(water_level - current_obj) * this->params.nlgd_params.nl_water_level_increase_pct;
        } else {
            water_level = water_level * std::exp(-this->params.nlgd_params.nl_water_level_decrease_exp * best_obj) + best_obj;
        }
    }

    template<typename Solution>
    bool NLGreatDeluge<Solution>::should_accept(double best_obj, double current_obj, double new_obj, double eps, Solution&, std::mt19937&) {
        this->current_obj = current_obj;
        this->best_obj = best_obj;
        return (new_obj < current_obj - eps || new_obj < water_level - eps);
    }
}

#endif