#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H

#include "../AcceptanceCriterion.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>

namespace mlpalns {
    /*! \brief This class models the Simulated Annealing acceptance criterion for ALNS */
    template<typename Solution>
    struct SimulatedAnnealing : public AcceptanceCriterion<Solution> {
        /*! Current temperature */
        double temperature;

        /*! Coefficient used to adjust the temperature in a geometric progression (i.e. exponential decay) */
        double alpha;

        /*! Step used to adjust the temperature in a linear progression (i.e. linear decrease) */
        double beta;

        /*! Start temperature */
        double start_temperature;

        /*! End temperature */
        double end_temperature;

        /*! Keep a local copy of best_obj to detect when it changes */
        double best_obj;

        /*! Keep a local copy of instance size */
        std::uint32_t instance_size;

        /*! Temperature last time we improved on the best solution */
        double last_improv_temperature;

        /*! Reheat every <m_iReheatingIterations> iterations */
        std::uint32_t reheating_iters;

        /*! Basic constructor */
        SimulatedAnnealing(Parameters& par) : AcceptanceCriterion<Solution>(par) {}

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
    void SimulatedAnnealing<Solution>::initialise() {
        // Read start and end temperature from the parameters and set current temp = start temp
        start_temperature = 0;
        end_temperature = 0;
        temperature = 1E-10;
        alpha = 1;
        beta = 0;
        best_obj = 0;
        last_improv_temperature = 0;
        reheating_iters = this->params.max_iters;
    }

    template<typename Solution>
    void SimulatedAnnealing<Solution>::calibrate(const CalibrationData& d) {
        instance_size = d.instance_size;

        if(this->params.sa_params.reheating_is_enabled) {
            // Temperature should go from min to max in this number of iterations
            reheating_iters = std::round(static_cast<long double>(this->params.max_iters) / static_cast<long double>(this->params.sa_params.reheating_times));

            // For weird values of iReheatingTimes, iReheatingIterations could
            // actaully be greater than the total number of remaining iterations!
            reheating_iters = std::min(reheating_iters, this->params.max_iters - this->params.prerun_iters);

            std::cout << "Reheating every " << reheating_iters << " iterations" << std::endl;
        } else {
            // No reheating is the same as reheating every... all the remaining iterations
            reheating_iters = this->params.max_iters - this->params.prerun_iters;
        }

        start_temperature = ((-this->params.sa_params.init_accept_ratio_50p * d.objective_value) / std::log(0.5)) /
                            std::pow(instance_size, this->params.sa_params.magic_number_exponent);
        end_temperature = ((-this->params.sa_params.end_accept_ratio_50p * d.objective_value) / std::log(0.5)) /
                          std::pow(instance_size, this->params.sa_params.magic_number_exponent);

        temperature = start_temperature;
        best_obj = d.objective_value;
        last_improv_temperature = temperature;

        if(end_temperature > start_temperature) {
            end_temperature = start_temperature;
        }

        // When decreasing temperature exponentially, use this multiplicative coefficient
        alpha = pow(end_temperature / start_temperature, 1.0 / reheating_iters);
        // When decreasing temeprature linearly, use this subtractive coefficient
        beta = (start_temperature - end_temperature) / reheating_iters;

        std::cout << "Start temp = " << start_temperature << ", end temp = " << end_temperature << std::endl;
        std::cout << "Initial exp coefficient: " << alpha << ", initial linear coefficient: " << beta << std::endl;
    }

    template<typename Solution>
    std::string SimulatedAnnealing<Solution>::get_print_str() const {
        return ("temp: " + std::to_string(static_cast<long double>(temperature)) + ", alpha: " + std::to_string(static_cast<long double>(alpha)) +
                ", beta: " + std::to_string(static_cast<long double>(beta)));
    }

    template<typename Solution>
    void SimulatedAnnealing<Solution>::update_parameters(std::uint32_t iter_number, double best_obj) {
        bool reheating = false;

        // Check if we have to reheat
        if(this->params.sa_params.reheating_is_enabled && (iter_number - this->params.prerun_iters) % reheating_iters == 0) {
            temperature = last_improv_temperature * this->params.sa_params.reheating_coefficient;
            reheating = true;
            std::cout << "Reheating! Iteration: " << iter_number << std::endl;
            std::cout << "New start temperature: " << temperature << std::endl;
        }

        if(this->params.sa_params.end_accept_ratio_refers_to_initial) {
            // End temperature depends on initial solution value =>
            // Adjust only when reheating

            if(reheating) {
                // Number of iterations in the current "sub-run" before reheating
                // This is equal to reheating_iters, unless there are fewer than
                // this number of iterations remaining, in which case we just use all
                // of them to go from start to end temperature.
                auto reheatingIterations = std::min(reheating_iters, this->params.max_iters - iter_number);

                std::cout << "Recalibrating parameters to go from " << temperature << " to " << end_temperature << " in " << reheatingIterations
                          << " iterations" << std::endl;

                alpha = std::min(pow(end_temperature / temperature, 1.0 / reheatingIterations), 1.0);
                beta = std::max((temperature - end_temperature) / reheatingIterations, 0.0);
            }
        } else {
            // End temperature depends on current best solution value =>
            // Adjust when reheating and when best solution value changes

            if(reheating || std::abs(best_obj - best_obj) > 1e-7) {
                auto iterationsToReachEndTemp = 0u;
                auto reheatingIterations = std::min(reheating_iters, this->params.max_iters - iter_number);

                if(reheating) {
                    iterationsToReachEndTemp = reheatingIterations;
                } else {
                    auto iterationsInThisSubRun = (iter_number - this->params.prerun_iters) % reheating_iters;

                    if(iter_number + reheatingIterations - iterationsInThisSubRun <= this->params.max_iters) {
                        iterationsToReachEndTemp = reheatingIterations - iterationsInThisSubRun;
                    } else {
                        iterationsToReachEndTemp = this->params.max_iters - iter_number;
                    }
                }

                best_obj = best_obj;
                end_temperature =
                    ((-this->params.sa_params.end_accept_ratio_50p * best_obj) / log(0.5)) / pow(instance_size, this->params.sa_params.magic_number_exponent);

                alpha = std::min(pow(end_temperature / temperature, 1.0 / iterationsToReachEndTemp), 1.0);
                beta = std::max((temperature - end_temperature) / iterationsToReachEndTemp, 0.0);
            }
        }

        // Update temperature
        if(!reheating) {
            if(this->params.sa_params.temperature_decrease_is_linear) {
                temperature -= beta;
            } else {
                temperature *= alpha;
            }
        }
    }

    template<typename Solution>
    bool SimulatedAnnealing<Solution>::should_accept(double best_obj, double current_obj, double new_obj, double eps, Solution&, std::mt19937& mt) {
        if(new_obj < best_obj) {
            last_improv_temperature = temperature;
        }

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        return (dist(mt) < exp((current_obj - new_obj) / temperature) - eps);
    }
}

#endif