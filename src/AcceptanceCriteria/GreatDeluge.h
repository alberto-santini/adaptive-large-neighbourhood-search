#ifndef GREAT_DELUGE_H
#define GREAT_DELUGE_H

#include "../AcceptanceCriterion.h"
#include <limits>

namespace mlpalns {
    /*! \brief This class models the Great Deluge acceptance criterion for ALNS */
    template<typename Solution>
    struct GreatDeluge : public AcceptanceCriterion<Solution> {
        /*! The water level. We accept only "underwater" solutions! */
        double water_level;

        /*! Objective value of the last current solution */
        double last_obj_value;

        /*! Basic constructor */
        GreatDeluge(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void GreatDeluge<Solution>::initialise() {
        water_level = std::numeric_limits<double>::max();
        last_obj_value = 0;
    }

    template<typename Solution>
    void GreatDeluge<Solution>::calibrate(const CalibrationData& d) {
        water_level = d.objective_value * this->params.gd_params.initial_water_level_ratio;
        last_obj_value = d.objective_value;
    }

    template<typename Solution>
    std::string GreatDeluge<Solution>::get_print_str() const {
        std::string wl = "";

        if(water_level > std::numeric_limits<double>::max() / 2) {
            wl = "infinite";
        } else {
            wl = std::to_string(static_cast<long double>(water_level));
        }
        return ("water level: " + wl);
    }

    template<typename Solution>
    void GreatDeluge<Solution>::update_parameters(std::uint32_t, double) {
        water_level -= this->params.gd_params.water_level_decrease_percentage * (water_level - last_obj_value);
    }

    template<typename Solution>
    bool GreatDeluge<Solution>::should_accept(double, double current_obj, double new_obj, double eps, Solution&, std::mt19937&) {
        last_obj_value = current_obj;
        return (new_obj < water_level - eps);
    }
}

#endif