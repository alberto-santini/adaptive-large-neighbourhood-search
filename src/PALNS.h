#ifndef __PALNS_H
#define __PALNS_H

#include "DestroyMethod.h"
#include "InitialSolutionCreator.h"
#include "Parameters.h"
#include "RepairMethod.h"
#include "AlgorithmVisitor.h"

#include "AcceptanceCriteria/ConservativeWorseAccept.h"
#include "AcceptanceCriteria/DiscreetWorseAccept.h"
#include "AcceptanceCriteria/GreatDeluge.h"
#include "AcceptanceCriteria/HillClimbing.h"
#include "AcceptanceCriteria/LateAcceptanceHillClimbing.h"
#include "AcceptanceCriteria/NLGreatDeluge.h"
#include "AcceptanceCriteria/RandomWalk.h"
#include "AcceptanceCriteria/RecordToRecordTravel.h"
#include "AcceptanceCriteria/SimulatedAnnealing.h"
#include "AcceptanceCriteria/ThresholdAcceptance.h"
#include "AcceptanceCriteria/WorseAccept.h"
#include "AcceptanceCriterion.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <chrono>

namespace mlpalns {
    /*! @brief This class models the PALNS algorithm solver
     *
     *  @tparam ProblemInstance Class which represents an instance of the problem.
     *                          It must implement the following methods:
     *                          * std::uint32_t getInstanceSize() giving a number representing the instance size.
     *  @tparam Solution        Class which represents a solution.
     *                          It must implement the following methods:
     *                          * Solution(const ProblemInstance&) an empty-solution constructor.
     *                          * Solution(const Solution&) a copy constructor. Can be the implicit default.
     *                          * void operator=(const Solution& other) an assignment operator. Can be the implicit default.
     *                          * double getCost() const, returning the cost of the solution (the smaller the better).
     */
    template<class ProblemInstance, class Solution>
    class PALNS {
    public:
        /*! @brief  Constructs the solver, starting from an instance of the problem.
         *          The initial, current, and best solutions are just default-constructed.
         *
         *  @param  instance The problem instance.
         */
        PALNS(const ProblemInstance& instance)
            : problem_instance(instance),
              algorithm_visitor(nullptr),
              acceptance_criterion(nullptr),
              best_solution(instance),
              current_solution(instance),
              initial_solution(instance) {}

        /*! @brief  Given a pointer to a destroy method and a string to describe it,
         *          it adds the method to the list of destroy methods.
         *
         *  @param   destroy_method    The destroy method to add.
         *  @param   description       A string description of the destroy method.
         *  @returns                   The index of the newly added method in the destroy methods vector.
         */
        std::size_t add_destroy_method(DestroyMethod<Solution>& destroy_method, std::string description) {
            destroy_methods.push_back(&destroy_method);
            destroy_methods_descriptions.push_back(description);

            return destroy_methods.size() - 1;
        }

        /*! @brief  Given a pointer to a repair method and a string to describe it,
         *          it adds the method to the list of repair methods. It is assumed that
         *          the repair method is compatible with all destroy methods. If this
         *          is not the case, use the overload of this method which takes three
         *          arguments.
         *
         *  @param   repair_method     The repair method to add.
         *  @param   description       A string description of the repair method.
         *  @returns                   The index of the newly added method in the repair methods vector.
         */
        std::size_t add_repair_method(RepairMethod<Solution>& repair_method, std::string description) {
            repair_methods.push_back(&repair_method);
            repair_methods_descriptions.push_back(description);
            repair_compatible_with_all.push_back(true);
            repair_compatibility_list.push_back(std::vector<std::size_t>());

            return repair_methods.size() - 1;
        }

        /*! @brief  Given a pointer to a repair method and a string to describe it,
         *          it adds the method to the list of repair methods. It is assumed that
         *          the repair method is only compatible with some of the destroy methods.
         *          The indices of the compatible destroy methods is passed in as the
         *          last argument.
         *
         *  @param  repair_method    The repair method to add.
         *  @param  description      A string description of the repair method.
         *  @param  compatible_with  A list of indices of destroy methods, with which the repair method is compatible.
         *  @returns                 The index of the newly added method in the repair methods vector.
         */
        std::size_t add_repair_method(RepairMethod<Solution>& repair_method, std::string description, std::vector<std::size_t> compatible_with) {
            repair_methods.push_back(&repair_method);
            repair_methods_descriptions.push_back(description);
            repair_compatible_with_all.push_back(false);
            repair_compatibility_list.push_back(compatible_with);

            return repair_methods.size() - 1;
        }

        /*! @brief  Sets the algorithm visitor, which is called at each iteration
         *          of the algorithm, and receives the algorithm status informations.
         *
         *  @param algorithm_visitor    The algorithm visitor.
         */
        void set_algorithm_visitor(std::unique_ptr<AlgorithmVisitor<Solution>>& algorithm_visitor) {
            this->algorithm_visitor = std::move(algorithm_visitor);
        }

        /*! @brief  Sets the acceptance criterion manually. For this acceptance
         *          criterion to be actually used, the user must set property
         *          \ref Parameters::acceptance_criterion_id to the value
         *          \ref Parameters::AcceptanceCriterionId::Custom.
         *
         *  @param acceptance_criterion The acceptance criterion to use.
         */
        void set_acceptance_criterion(std::unique_ptr<AcceptanceCriterion<Solution>>& acceptance_criterion) {
            this->acceptance_criterion = std::move(acceptance_criterion);
        }

        /*! @brief  Starts the solution algorithm.
         *
         *  @param   start_solution    The initial solution.
         *  @param   num_threads       Number of parallel threads to use.
         *  @param   parameters        Algorithm parameters.
         *  @returns                   The best solution encountered.
         */
        Solution go(const Solution& start_solution, std::uint32_t num_threads, Parameters parameters = Parameters()) {
            // Clear up all the local data
            reset_local_parameters(start_solution, parameters);

            // Initialise method-specific parameters in the acceptance criterion
            acceptance_criterion->initialise();

            // Invoke the algorithm visitor at the beginning of the algorithm
            algorithm_visitor->on_algorithm_start(
                destroy_methods,
                repair_methods,
                destroy_methods_descriptions,
                repair_methods_descriptions
            );

            // Make a preliminary run
            do_preliminary_run(num_threads);

            // Record when we started
            start_time = std::chrono::high_resolution_clock::now();

            // This is the real computation:
            start_threads(num_threads);

            return best_solution;
        }

        /*! @brief  Repeats the solution algorithm many times, e.g. for statistical testing purposes.
         *
         *  @param   start_solution    The initial solution
         *  @param   retries           Number of times the algorithm should be executed
         *  @param   num_threads       Number of parallel threads to use, for each run of the algorithm
         *  @param   parameters        Algorithm parameters
         *  @returns                   The best solution encountered over all the reruns
         */
        Solution repeat_solver(const Solution& start_solution, std::uint32_t retries, std::uint32_t num_threads, Parameters parameters = Parameters()) {
            Solution best_solution(start_solution);

            // Repeat the solution process for retries times
            for(auto i = 0u; i < retries; ++i) {
                Solution solution(start_solution);

                // Launch the PALNS algorithm
                go(solution, num_threads, parameters);

                // Update the best solution of all re-runs
                if(solution.getCost() < best_solution.getCost()) {
                    best_solution = solution;
                }
            }

            return best_solution;
        }

        /*! @brief  Repeats the solution algorithm many times, e.g. for statistical testing purposes.
         *
         *  @param initial_solution_generator    A pointer to a solution generator, used to obtain the initial solution.
         *  @param retries                       Number of times the algorithm should be executed.
         *  @param num_threads                   Number of parallel threads to use, for each run of the algorithm.
         *  @param parameters                    Algorithm parameters.
         *  @returns                             The best solution encountered over all the reruns.
         */
        Solution repeat_solver(InitialSolutionCreator<Solution, ProblemInstance>* initial_solution_generator, std::uint32_t retries, std::uint32_t num_threads, Parameters parameters = Parameters()) {
            std::mt19937::result_type random_data[std::mt19937::state_size];
            std::random_device source;
            std::generate(std::begin(random_data), std::end(random_data), std::ref(source));
            std::seed_seq seeds(std::begin(random_data), std::end(random_data));
            std::mt19937 mt(seeds);

            return repeat_solver(initial_solution_generator->create_initial_solution(problem_instance, mt), retries, num_threads, parameters);
        }

    private:
        /*! @brief  Fire up a single thread.
         *
         *  @param   thread_id     Progressive id of the thread.
         */
        void start_thread(std::uint32_t thread_id) {
            if(!acceptance_criterion) {
                throw std::runtime_error("Cannot start PALNS without an acceptance criterion!");
            }

            // Generate random seed for this thread
            std::mt19937::result_type random_data[std::mt19937::state_size];
            std::iota(std::begin(random_data), std::end(random_data), thread_id);
            std::seed_seq seeds(std::begin(random_data), std::end(random_data));
            std::mt19937 mt = std::mt19937(seeds);

            // Local copies of destroy and repair methods
            std::vector<std::unique_ptr<DestroyMethod<Solution>>> loc_destroy_methods;
            std::vector<std::unique_ptr<RepairMethod<Solution>>> loc_repair_methods;
            std::vector<std::vector<std::size_t>> loc_repair_compatibility;

            auto make_local_copies_of_dr_methods = [&,this] () -> void {
                loc_destroy_methods.clear();
                loc_repair_methods.clear();

                std::lock_guard<std::mutex> _(destroy_repair_mtx);

                for(auto& method : destroy_methods) {
                    loc_destroy_methods.push_back(method->clone());
                }

                loc_repair_compatibility.resize(destroy_methods.size());
                for(auto i = 0u; i < repair_methods.size(); ++i) {
                    loc_repair_methods.push_back(repair_methods[i]->clone());

                    if(!repair_compatible_with_all[i]) {
                        const auto& comp_vec = repair_compatibility_list[i];
                        for(auto j = 0u; j < comp_vec.size(); ++j) {
                            loc_repair_compatibility[comp_vec[j]].push_back(i);
                        }
                    } else {
                        for(auto dest_id = 0u; dest_id < destroy_methods.size(); ++dest_id) {
                            loc_repair_compatibility[dest_id].push_back(i);
                        }
                    }
                }
            };

            make_local_copies_of_dr_methods();

            // Counter of number of iterations the current thread has worked since last time it synchronised its data
            std::uint64_t local_iterations_count = 0;

            // Last iteration in which we recorded that the number of iterations without improvement was too high,
            // and we called the algorithm visitor to do something about it.
            std::uint32_t last_iter_global_alarm_visitor_called = 0u;

            // Create a temporary solution to apply destroy and repair to
            Solution tmp_sol = threadsafe_clone_solution(current_solution, current_solution_mtx);

            // Create a local solution (will be updated if tmp_sol is accepted)
            Solution loc_current_solution(tmp_sol);

            // ...and get its objective value
            double local_best_obj = tmp_sol.getCost();

            // This is true when we finish! :)
            bool done = false;

            std::uniform_real_distribution<float> dist(0.0, 1.0);

            // Main loop
            while(!done) {
                // ****** 1. Check if we are done ******
                
                // Local current iteration
                auto cur_iter = 0u;

                {
                    // Acquire lock on iteration counter
                    std::lock_guard<std::mutex> _(num_iter_mtx);

                    // We increase iteration counter by one in advance so that we are sure that do not perform more than max_iters iterations.
                    if(num_iter < max_iters) {
                        num_iter++;
                    } else {
                        done = true;
                    }

                    // Check for timeout
                    auto current_time = std::chrono::high_resolution_clock::now();
                    auto elapsed_time = std::chrono::duration_cast<std::chrono::duration<double>>(current_time - start_time).count();
                    if(elapsed_time > params.max_seconds) {
                        done = true;
                    }

                    cur_iter = num_iter;
                }

                // ****** 2. If we are not done, let's rock! ******
                if(!done) {
                    // Check if we should return the best known solution
                    if(dist(mt) < params.return_to_best_known_prob) {
                        {
                            std::lock_guard<std::mutex> _(best_solution_mtx);
                            tmp_sol = best_solution;
                        }
                    }

                    double best_sol_cost = threadsafe_get_solution_cost(best_solution, best_solution_mtx);

                    // Check if it's time to print output
                    if((cur_iter % print_output_iter) == 0) {
                        std::string acceptance_data;

                        // Get data from the acceptance criterion
                        {
                            std::lock_guard<std::mutex> _(acceptance_criterion_mtx);
                            acceptance_data = acceptance_criterion->get_print_str();
                        }

                        std::cout << cur_iter << ", thread " << thread_id << ", best solution: " << best_sol_cost << ", current solution: " << tmp_sol.getCost()
                                  << ", " << acceptance_data << std::endl;
                    }

                    std::size_t chosen_destroy_id = 0, chosen_repair_id = 0;

                    // Select destroy and repair heuristics
                    {
                        std::lock_guard<std::mutex> _(weights_mtx);
                        chosen_destroy_id = roulette_wheel_selection(destroy_weights, mt);
                        chosen_repair_id = roulette_wheel_selection(loc_repair_compatibility[chosen_destroy_id], repair_weights, mt);
                    }

                    // Destroy solution
                    loc_destroy_methods[chosen_destroy_id]->destroy_solution(tmp_sol, mt);
                    // Repair solution
                    loc_repair_methods[chosen_repair_id]->repair_solution(tmp_sol, mt);

                    // Get new solution cost
                    double new_sol_cost = tmp_sol.getCost();

                    // Update acceptance criterion parameters
                    {
                        std::lock_guard<std::mutex> _(acceptance_criterion_mtx);
                        acceptance_criterion->update_parameters(cur_iter, best_sol_cost);
                    }

                    bool improved = false;
                    bool accepted = false;
                    double cur_solution_cost;

                    if(local_iterations_count < max_iters_before_sync) {
                        // If we should only update local data...

                        // Increment local iteration counter
                        local_iterations_count++;

                        // Get current objective solution
                        cur_solution_cost = loc_current_solution.getCost();

                        // Determine wether to accept the solution or not
                        {
                            std::lock_guard<std::mutex> _(acceptance_criterion_mtx);
                            accepted = acceptance_criterion->should_accept(best_sol_cost, cur_solution_cost, new_sol_cost, eps, tmp_sol, mt);
                        }

                        if(accepted) {
                            // Update current (local) solution
                            loc_current_solution = tmp_sol;

                            // Mark as improved, if we improved on the current solution
                            if(new_sol_cost < cur_solution_cost - eps) {
                                improved = true;
                            }

                            // Update the current objective solution
                            cur_solution_cost = loc_current_solution.getCost();
                        } else {
                            // Revert back tmp_sol to the current (local) solution
                            tmp_sol = loc_current_solution;
                        }
                    } else {
                        // If we should update shared data...

                        // Reset local iteration counter
                        local_iterations_count = 0;

                        {
                            // Acquire mutex on the current solution
                            std::lock_guard<std::mutex> _c(current_solution_mtx);

                            // Update the current objective value
                            cur_solution_cost = current_solution.getCost();

                            // Should local current solution be the global current solution?
                            // Only necessary to check if we use "local iterations"
                            // That is, if max_iters_before_sync > 0
                            // Check current local solution to be sure that it has a chance of getting accepted, even if
                            // The last solution found by the thread (new_sol_cost, see below) is maybe not so good
                            if(max_iters_before_sync > 0) {
                                bool tmp_accepted = false;

                                {
                                    std::lock_guard<std::mutex> _(acceptance_criterion_mtx);
                                    tmp_accepted = acceptance_criterion->should_accept(best_sol_cost, cur_solution_cost, loc_current_solution.getCost(), eps,
                                                                                       loc_current_solution, mt);
                                }

                                if(tmp_accepted) {
                                    current_solution = loc_current_solution;
                                    cur_solution_cost = current_solution.getCost();
                                }
                            }

                            // Check if the last solution found by the thread should be accepted
                            {
                                std::lock_guard<std::mutex> _(acceptance_criterion_mtx);
                                accepted = acceptance_criterion->should_accept(best_sol_cost, cur_solution_cost, new_sol_cost, eps, tmp_sol, mt);
                            }

                            // Determine wether to accept the solution or not
                            if(accepted) {
                                // Update current (shared) solution
                                current_solution = tmp_sol;

                                // Mark as improved, if we improved on the current solution
                                if(new_sol_cost < cur_solution_cost - eps) {
                                    improved = true;
                                }

                                // Update the current objective solution
                                cur_solution_cost = current_solution.getCost();
                            } else {
                                // Revert back tmp_sol to the current (shared) solution
                                tmp_sol = current_solution;
                            }

                            // Are we using the local current solution? If so, it has to be updated
                            if(max_iters_before_sync > 0) {
                                loc_current_solution = current_solution;
                            }
                        } // Release mutex on the current solution
                    }

                    // Let's check if we found the new global best
                    bool new_global_best = false;

                    if(accepted) {
                        // Check if best solution should be updated:
                        // First check the local best solution, if the current solution is better than that, then also check the
                        // global best one. We really didn't need the local best solution, but having it saves us from some lock
                        // operations (quite a few I would think).

                        // If new obj improves on local best
                        if(new_sol_cost < local_best_obj - eps) {
                            // Update local best
                            local_best_obj = new_sol_cost;

                            {
                                std::lock_guard<std::mutex> _(best_solution_mtx);
                                if(new_sol_cost < best_solution.getCost() - eps) {
                                    // Updated global best
                                    best_solution = tmp_sol;
                                    last_gobal_improvement_iter = num_iter;
                                    new_global_best = true;
                                }
                            }
                        }

                        // If we actually did better than the global best, we print it
                        if(new_global_best) {
                            std::cout.precision(10);
                            std::cout << cur_iter << ", thread " << thread_id << ", New best solution: " << tmp_sol.getCost() << ". ";
                            std::cout << "Destroy method: " << destroy_methods_descriptions[chosen_destroy_id]
                                      << ". Repair method: " << repair_methods_descriptions[chosen_repair_id] << std::endl;
                        }
                    }

                    // Calculate new score
                    float new_score = calculate_score(accepted, improved, new_global_best);

                    // Update weights
                    {
                        std::lock_guard<std::mutex> _(weights_mtx);
                        destroy_weights[chosen_destroy_id] =
                            destroy_weights[chosen_destroy_id] * params.score_decay +
                            new_score * (1 - params.score_decay);
                        repair_weights[chosen_repair_id] =
                            repair_weights[chosen_repair_id] * params.score_decay +
                            new_score * (1 - params.score_decay);
                    }

                    // Do we have an algorithm visitor? If so, call it
                    if(algorithm_visitor && num_iter > params.prerun_iters) {
                        AlgorithmStatus<Solution> status(params, tmp_sol, best_solution);
                        bool should_raise_global_improvement_alarm = false;

                        {
                            std::lock_guard<std::mutex> _ni(num_iter_mtx);
                            status.iter_number = num_iter - params.prerun_iters;

                            // If the last time the global best was improved happend iters_without_improvement_alarm iterations ago:
                            if(num_iter >= last_gobal_improvement_iter + params.iters_without_improvement_alarm) {
                                // If we haven't raised this alarm during the past iters_without_improvement_alarm previous iterations:
                                if(last_iter_global_alarm_visitor_called <= num_iter - params.iters_without_improvement_alarm) {
                                    // We should raise the alarm!
                                    should_raise_global_improvement_alarm = true;
                                    last_iter_global_alarm_visitor_called = num_iter;
                                }
                            }
                        }

                        status.destroy_method_id = chosen_destroy_id;
                        status.repair_method_id = chosen_repair_id;
                        status.accepted = accepted;
                        status.new_best = new_global_best;

                        {
                            // The visitor could modify the best solution:
                            std::lock_guard<std::mutex> _bs(best_solution_mtx);
                            algorithm_visitor->on_iteration_end(status);
                        }

                        if(should_raise_global_improvement_alarm) {
                            // The visitor could modify the detroy/repair methods:
                            {
                                std::lock_guard<std::mutex> _dr(destroy_repair_mtx);
                                algorithm_visitor->on_many_iters_without_improvement(destroy_methods, repair_methods);
                            }
                            make_local_copies_of_dr_methods();
                        }
                    }
                }
            }
        }

        /*! @brief  Fire up many threads.
         *
         *  @param   num_threads   Number of threads to start
         */
        void start_threads(std::uint32_t num_threads) {
            std::vector<std::thread> threads(num_threads);
            for(auto i = 0u; i < num_threads; ++i) {
                threads[i] = std::thread([i, this]() { start_thread(i); });
            }

            for(auto& thread : threads) {
                thread.join();
            }
        }

        /*! @brief  Calculates the score of a destroy/repair method which produced a new solution
         *          that was accepted / improved on current / improved on best.
         *
         *  @param   is_accepted           The produced solution was accepted?
         *  @param   is_improved           The produced solution improved on the current?
         *  @param   is_new_global_best    The produced solution improved on the best?
         *  @return                        The new score.
         */
        float calculate_score(bool is_accepted, bool is_improved, bool is_new_global_best) const {
            float score = 1.0f;

            // Increase score if accepted
            if(is_accepted) {
                score = std::max(score, params.score_mult_accepted);
            }

            // Increase score if improved on current
            if(is_improved) {
                score = std::max(score, params.score_mult_improved);
            }

            // Increase score if improved on overall best
            if(is_new_global_best) {
                score = std::max(score, params.score_mult_global_best);
            }

            return score;
        }

        /*! @brief  Performs a roulette-wheel selection which returns a random index, with
         *          probability directly proportional to its weight.
         *
         *  @param  weights   Vector of weights.
         *  @param  mt        Random numbers generator.
         *  @return           The index of the chosen element.
         */
        std::size_t roulette_wheel_selection(const std::vector<float>& weights, std::mt19937& mt) const {
            // Sum of all weights
            float sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0f);
            std::uniform_real_distribution<float> dist(0.0, sum_weights);

            // Pick a random number between 0 and the sum of all weights
            float rnd = dist(mt);

            // Now look for the right weight
            float summed_weights = 0;
            std::size_t method = weights.size() - 1; // This will store the chosen index

            for(std::size_t i = 0; i < weights.size() - 1; ++i) {
                summed_weights += weights[i];

                if(rnd <= summed_weights) {
                    method = i;
                    break;
                }
            }

            return method;
        }

        /*! @brief  Performs a roulette-wheel selection which returns a random index, with
         *          probability directly proportional to its weight.
         *
         *  @param  compatible_methods    Vector of indices to chose from.
         *  @param  weights               Vector of weights.
         *  @param  mt                    Random numbers generator.
         *  @return                       The index of the chosen element.
         */
        std::size_t roulette_wheel_selection(const std::vector<std::size_t>& compatible_methods, const std::vector<float>& weights, std::mt19937& mt) const {
            // Like the previous method, but here only certain elements can be chosen, because they correspond to repair methods that are compatible with the chosen
            // destroy method

            assert(!compatible_methods.empty());

            float sum_weights = 0.0;

            for(auto method_id : compatible_methods) {
                sum_weights += weights[method_id];
            }

            std::uniform_real_distribution<float> dist(0.0, sum_weights);
            float rnd = dist(mt);

            float summed_weights = 0;
            std::size_t method = compatible_methods.back();

            for(auto method_id : compatible_methods) {
                summed_weights += weights[method_id];

                if(rnd <= summed_weights) {
                    method = method_id;
                    break;
                }
            }

            return method;
        }

        /*! @brief  Brings back the internal state of PALNS to its starting point.
         *
         *  @param  start_solution    The initial solution.
         *  @param  parameters        Algorithm parameters.
         */
        void reset_local_parameters(Solution start_solution, Parameters& parameters)  {
            // We start at iteration 0
            num_iter = 0u;

            // No improvement yet.
            last_gobal_improvement_iter = 0u;

            // Reset the start time
            start_time = std::chrono::high_resolution_clock::now();

            // Save the starting solution
            this->initial_solution = start_solution;

            // Best solution is the starting solution
            best_solution = start_solution;

            // Current solution is the starting solution
            current_solution = start_solution;

            // Import parameters for local use
            max_iters_before_sync = parameters.max_independent_thread_iters;
            max_iters = parameters.max_iters;

            // Save local copy of parameters file
            params = parameters;

            // Shared parameters for adaptive weight adjustment
            destroy_weights.clear();
            destroy_weights.resize(destroy_methods.size(), 1.0f);
            repair_weights.clear();
            repair_weights.resize(repair_methods.size(), 1.0f);

            // Create a new acceptance criterion object
            switch(params.acceptance_criterion_id) {
                case Parameters::AcceptanceCriterionId::HillClimbing:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new HillClimbing<Solution>(params));
                    std::cout << "Acceptance criterion: Hill climbing" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::SimulatedAnnealing:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new SimulatedAnnealing<Solution>(params));
                    std::cout << "Acceptance criterion: Simulated Annealing" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::ThresholdAcceptance:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new ThresholdAcceptance<Solution>(params));
                    std::cout << "Acceptance criterion: Threshold acceptance" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::GreatDeluge:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new GreatDeluge<Solution>(params));
                    std::cout << "Acceptance criterion: Great Deluge" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::RecordToRecordTravel:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new RecordToRecordTravel<Solution>(params));
                    std::cout << "Acceptance criterion: Record-To-Record Travel" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::LateAcceptanceHillClimbing:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new LateAcceptanceHillClimbing<Solution>(params));
                    std::cout << "Acceptance criterion: Late Acceptance Hill Climbing" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::NonLinearGreatDeluge:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new NLGreatDeluge<Solution>(params));
                    std::cout << "Acceptance criterion: Non-linear Great Deluge" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::WorseAccept:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new WorseAccept<Solution>(params));
                    std::cout << "Acceptance criterion: Worse-accept" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::ConservativeWorseAccept:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new ConservativeWorseAccept<Solution>(params));
                    std::cout << "Acceptance criterion: Conservative Worse-accept" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::DiscreetWorseAccept:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new DiscreetWorseAccept<Solution>(params));
                    std::cout << "Acceptance criterion: Discreet Worse-accept" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::RandomWalk:
                    acceptance_criterion = std::unique_ptr<AcceptanceCriterion<Solution>>(new RandomWalk<Solution>(params));
                    std::cout << "Acceptance criterion: Random walk" << std::endl;
                    break;

                case Parameters::AcceptanceCriterionId::Custom:
                    std::cout << "Acceptance criterion: Custom" << std::endl;
                    break;
            }
        }

        /*! @brief  Executes a preliminary run, used for calibration of the acceptance criteria.
         *
         *  @param  num_threads   Number of parallel threads to use.
         */
        void do_preliminary_run(std::uint32_t num_threads) {
            // How many iterations in the preliminary run?
            max_iters = params.prerun_iters;

            if(max_iters > 0u) {
                // Save the original number of iterations without synchronising
                auto max_iters_before_sync_bak = max_iters_before_sync;

                // Set to 0 the number of iterations without synchronising
                max_iters_before_sync = 0u;

                // Tell the acceptance criterion that we are running a preliminary run
                acceptance_criterion->set_preliminary_run();

                // Perform the preliminary run
                start_threads(num_threads);

                // Tell the acceptance criterion the preliminary run has finished
                acceptance_criterion->unset_preliminary_run();

                // Restore the original number of iterations without synchronising
                max_iters_before_sync = max_iters_before_sync_bak;
            }

            // Get the current solution's objective value
            double current_obj = best_solution.getCost();

            // Set the current solution
            current_solution = best_solution;

            // Get the instance size
            std::uint32_t instance_size = static_cast<std::uint32_t>(problem_instance.getInstanceSize());

            // How many iterations to perform in the normal run?
            max_iters = params.max_iters;

            // Create the data to be used to calibrate the acceptance criterion
            CalibrationData calibration_data(current_obj, instance_size, max_iters);

            // Perform the calibration!
            acceptance_criterion->calibrate(calibration_data);

            // Invoke the appropriate algorithm visitor
            algorithm_visitor->on_prerun_end(destroy_methods, repair_methods);
        }

        /*! @brief Clones a solution (but first locks the corresponding mutex).
         *
         *  @param sol       The solution to clone.
         *  @param sol_mutex The mutex protecting the solution.
         *  @returns         The new, cloned solution.
         */
        Solution threadsafe_clone_solution(const Solution& sol, std::mutex& sol_mutex) {
            std::lock_guard<std::mutex> _(sol_mutex);
            return Solution(sol);
        }

        /*! @brief Gets a solution's cost (but first locks the corresponding mutex).
         *
         *  @param sol       The solution of which we need the cost.
         *  @param sol_mutex The mutex protecting the solution.
         *  @return          The cost of the solution.
         */
        double threadsafe_get_solution_cost(const Solution& sol, std::mutex& sol_mutex) {
            std::lock_guard<std::mutex> _(sol_mutex);
            return sol.getCost();
        }

        /*! @brief Instance of the problem. */
        const ProblemInstance problem_instance;

        /*! @brief Names of destroy heuristics. */
        std::vector<std::string> destroy_methods_descriptions;

        /*! @brief Names of repair heuristics. */
        std::vector<std::string> repair_methods_descriptions;

        /*! @brief  List of repair methods which are guaranteed to be compatible
         *          with any destroy method. */
        std::deque<bool> repair_compatible_with_all;

        /*! @brief Compatibility table for repair and destroy methods.
         *
         *         The vector at index [destroy] contains a list of indices of
         *         compatible repair methods.
         */
        std::vector<std::vector<std::size_t>> repair_compatibility_list;

        // ----- Member variables that are copied in each thread. -----

        /*! @brief List of destroy methods. */
        std::vector<DestroyMethod<Solution>*> destroy_methods;

        /*! @brief List of repair methods. */
        std::vector<RepairMethod<Solution>*> repair_methods;

        // ----- Shared variables. Each is protected by a mutex. -----

        /*! @brief Best solution encountered so far. */
        Solution best_solution;

        /*! @brief Current solution. */
        Solution current_solution;

        /*! @brief Initial solution. */
        Solution initial_solution;

        /*! @brief Current iteration number. */
        std::uint32_t num_iter;

        /*! @brief Last iteration at which the global best was improved. */
        std::uint32_t last_gobal_improvement_iter;

        /*! @brief Start time of the algorithm. */
        std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

        /*! @brief Weights associated with the destroy methods. */
        std::vector<float> destroy_weights;

        /*! @brief Weights associated with the repair methods. */
        std::vector<float> repair_weights;

        /*! @brief Acceptance criterion to be used. */
        std::unique_ptr<AcceptanceCriterion<Solution>> acceptance_criterion;

        /*! @brief  Algorithm visitor to be used. Called at each iteration.
         *          Can be nullptr. */
        std::unique_ptr<AlgorithmVisitor<Solution>> algorithm_visitor;

        // ----- Mutexes -----

        /*! @brief Mutex on the current solution. */
        std::mutex current_solution_mtx;

        /*! @brief Mutex on the best solution (so far). */
        std::mutex best_solution_mtx;

        /*! @brief Mutex on the iterations counter. */
        std::mutex num_iter_mtx;

        /*! @brief Mutex on the acceptance criterion. */
        std::mutex acceptance_criterion_mtx;

        /*! @brief Mutex on the weights of destroy/repair methods. */
        std::mutex weights_mtx;

        /*! @brief Mutex on the destroy/repair methods. */
        std::mutex destroy_repair_mtx;

        // ----- Shared constants. These must not be changed while PALNS is running. -----

        /*! @brief Algorithm parameters. */
        Parameters params;

        /*! @brief Maximum iterations to perform. */
        std::uint32_t max_iters;

        /*! @brief  Number of iterations each thread should advance before syncing
         *          with the other threads. */
        std::uint32_t max_iters_before_sync;

        // ----- Static constexpr members ------

        /*! @brief Precision to use when comparing solution scores. */
        static constexpr double eps = 1e-6;

        /*! @brief  We print out updates on the solution process every 1000
         *          iterations. */
        static constexpr std::uint32_t print_output_iter = 1000;
    };
}

#endif
