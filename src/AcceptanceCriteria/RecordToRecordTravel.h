#ifndef RECORD_DEVIATION_H
#define RECORD_DEVIATION_H

#include "../AcceptanceCriterion.h"

namespace mlpalns {
    /*! \brief This class models the Record Deviation acceptance criterion for ALNS */
    template<typename Solution>
    struct RecordToRecordTravel : public AcceptanceCriterion<Solution> {
        /*! Deviation */
        double deviation;

        /*! Deviation step, when decrease is linear */
        double deviation_step;

        /*! Deviation multiplicative coefficient, when decrease is exponential */
        double deviation_exp_coeff;

        /*! Basic constructor */
        RecordToRecordTravel(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void RecordToRecordTravel<Solution>::initialise() {
        deviation = 1.0;
        deviation_step = 0.0;
        deviation_exp_coeff = 1.0;
    }

    template<typename Solution>
    void RecordToRecordTravel<Solution>::calibrate(const CalibrationData&) {
        deviation = this->params.rrt_params.start_deviation;
        deviation_step =
            (this->params.rrt_params.start_deviation - this->params.rrt_params.end_deviation) / (this->params.max_iters - this->params.prerun_iters);
        deviation_exp_coeff =
            pow(this->params.rrt_params.start_deviation / this->params.rrt_params.end_deviation, 1.0 / (this->params.max_iters - this->params.prerun_iters));
    }

    template<typename Solution>
    std::string RecordToRecordTravel<Solution>::get_print_str() const {
        return ("deviation: " + std::to_string(static_cast<long double>(deviation)));
    }

    template<typename Solution>
    void RecordToRecordTravel<Solution>::update_parameters(std::uint32_t, double) {
        if(this->params.rrt_params.deviation_decrease_is_linear) {
            deviation -= deviation_step;
        } else {
            deviation *= deviation_exp_coeff;
        }
    }

    template<typename Solution>
    bool RecordToRecordTravel<Solution>::should_accept(double best_obj, double, double new_obj, double eps, Solution&, std::mt19937&) {
        double gap = (new_obj - best_obj) / new_obj;
        return (gap < deviation - eps);
    }
}

#endif