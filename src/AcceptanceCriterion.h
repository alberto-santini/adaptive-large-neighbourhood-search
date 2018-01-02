#ifndef ACCEPTANCE_CRITERION_H
#define ACCEPTANCE_CRITERION_H

#include "Parameters.h"

#include <random>
#include <string>

namespace mlpalns {
    /*! \brief This class encapsulate the data that the acceptance criterion needs to auto-calibrate its parameters */
    struct CalibrationData {
        /*! The value of a solution */
        double objective_value;

        /*! The size of the instance */
        std::uint32_t instance_size;

        /*! Maximum number of iterations (if this is the stopping criterion) */
        std::uint32_t max_iter;

        /*! Empty constructor */
        CalibrationData() {}

        /*! Basic constructor */
        CalibrationData(double obj, std::uint32_t si, std::uint32_t mi) : objective_value(obj), instance_size(si), max_iter(mi) {}
    };

    /*! \brief This class models a base, generic acceptance criterion */
    template<typename Solution>
    struct AcceptanceCriterion {
        /*! Reference to parameters object */
        Parameters& params;

        /*! Stateful parameter that is true iff the acceptance criterion is being used in a ``preliminary run'' */
        bool preliminary_run;

        /*! Basic constructor */
        AcceptanceCriterion(Parameters& par) : params(par), preliminary_run(false) {}

        /*! Setter for the preliminary run flag */
        void set_preliminary_run() { preliminary_run = true; }

        /*! Unsetter for the preliminary run flag */
        void unset_preliminary_run() { preliminary_run = false; }

        // Methods each subclass should implement:

        /*! This method initialise method-specific parameters */
        virtual void initialise() = 0;

        /*! This method is used to calibrate the parameters, given a starting solution's data
         *
         *  @param d is the data to be used for calibration
         */
        virtual void calibrate(const CalibrationData& d) = 0;

        /*! This method provides some info about the current status of the acceptance criterion, to be printed out during the computation */
        virtual std::string get_print_str() const = 0;

        /*! This method updates the acceptance criterion's parameters, based on running info
         *
         *  @param iter_number is the current iteration number
         *  @param best_obj is the value of the current best solution
         */
        virtual void update_parameters(std::uint32_t iter_number, double best_obj) = 0;

        /*! This method returns true iff the solution should be accepted according to the acceptance criterion
         *
         *  @param best_obj is the objective value of the best solution (globally) found so far
         *  @param current_obj is the objective value of the current solution
         *  @param new_obj is the objective value of the new solution (which we must determine wether to accept or not)
         *  @param eps is an epsilon to be used in numerical comaprisons
         *  @param new_sol is a reference to the new solution found
         *  @param mt is a reference to a random number generator
         */
        virtual bool should_accept(double best_obj, double current_obj, double new_obj, double eps, Solution& new_sol, std::mt19937& mt) = 0;
    };
}

#endif