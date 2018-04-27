#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdint>
#include <string>

#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>

namespace mlpalns {
    struct Parameters {
        struct SimulatedAnnealingParameters {
            double init_accept_ratio_50p;
            double end_accept_ratio_50p;
            double magic_number_exponent;
            bool end_accept_ratio_refers_to_initial;
            bool temperature_decrease_is_linear;
            bool reheating_is_enabled;
            double reheating_coefficient;
            double reheating_times;

            SimulatedAnnealingParameters() = default;

            SimulatedAnnealingParameters(double ia, double ea, double mne, bool earti, bool tdl, bool er, double rc, double rt)
                : init_accept_ratio_50p(ia),
                  end_accept_ratio_50p(ea),
                  magic_number_exponent(mne),
                  end_accept_ratio_refers_to_initial(earti),
                  temperature_decrease_is_linear(tdl),
                  reheating_is_enabled(er),
                  reheating_coefficient(rc),
                  reheating_times(rt) {}
        };

        struct ThresholdParameters {
            double start_threshold;
            double end_threshold;
            bool threshold_decrease_is_linear;

            ThresholdParameters() = default;

            ThresholdParameters(double st, double et, bool tdl) : start_threshold(st), end_threshold(et), threshold_decrease_is_linear(tdl) {}
        };

        struct GreatDelugeParameters {
            double initial_water_level_ratio;
            double water_level_decrease_percentage;

            GreatDelugeParameters() = default;

            GreatDelugeParameters(double wl, double inc) : initial_water_level_ratio(wl), water_level_decrease_percentage(inc) {}
        };

        struct RecordDeviationParameters {
            double start_deviation;
            double end_deviation;
            bool deviation_decrease_is_linear;

            RecordDeviationParameters() = default;

            RecordDeviationParameters(double sd, double ed, bool ddl) : start_deviation(sd), end_deviation(ed), deviation_decrease_is_linear(ddl) {}
        };

        struct LateAcceptanceParameters {
            std::uint32_t list_size;
            bool accept_non_worsening;

            LateAcceptanceParameters() = default;

            LateAcceptanceParameters(std::uint32_t ls, bool anw) : list_size(ls), accept_non_worsening(anw) {}
        };

        struct NLGreatDelugeParameters {
            double nl_initial_water_level_ratio;
            double nl_gap_to_increase_water_level;
            double nl_water_level_increase_pct;
            double nl_water_level_decrease_exp;

            NLGreatDelugeParameters() = default;

            NLGreatDelugeParameters(double iwl, double gti, double wlip, double wlde)
                : nl_initial_water_level_ratio(iwl),
                  nl_gap_to_increase_water_level(gti),
                  nl_water_level_increase_pct(wlip),
                  nl_water_level_decrease_exp(wlde) {}
        };

        struct WorseAcceptParameters {
            double start_prob;
            double end_prob;
            bool prob_decrease_is_linear;

            WorseAcceptParameters() = default;

            WorseAcceptParameters(double sp, double ep, bool lin) : start_prob(sp), end_prob(ep), prob_decrease_is_linear(lin) {}
        };

        struct ConservativeWorseAcceptParameters {
            double start_prob;
            double end_prob;
            bool prob_decrease_is_linear;

            ConservativeWorseAcceptParameters() = default;

            ConservativeWorseAcceptParameters(double sp, double ep, bool lin) : start_prob(sp), end_prob(ep), prob_decrease_is_linear(lin) {}
        };

        struct DiscreetWorseAcceptParameters {
            std::uint32_t consecutive_rejects_for_100p;
            bool prob_increase_is_linear;

            DiscreetWorseAcceptParameters() = default;

            DiscreetWorseAcceptParameters(std::uint32_t cr, bool lin) : consecutive_rejects_for_100p(cr), prob_increase_is_linear(lin) {}
        };

        enum class AcceptanceCriterionId : std::uint32_t {
            HillClimbing,
            SimulatedAnnealing,
            ThresholdAcceptance,
            GreatDeluge,
            RecordToRecordTravel,
            LateAcceptanceHillClimbing,
            NonLinearGreatDeluge,
            WorseAccept,
            ConservativeWorseAccept,
            DiscreetWorseAccept,
            RandomWalk,
            Custom
        };

        enum class AcceptanceParamsBase : std::uint32_t {
            Time, Iterations
        };

        AcceptanceCriterionId acceptance_criterion_id;
        AcceptanceParamsBase acceptance_params_base;

        SimulatedAnnealingParameters sa_params;
        ThresholdParameters ta_params;
        GreatDelugeParameters gd_params;
        RecordDeviationParameters rrt_params;
        LateAcceptanceParameters lahc_params;
        NLGreatDelugeParameters nlgd_params;
        WorseAcceptParameters wa_params;
        ConservativeWorseAcceptParameters cwa_params;
        DiscreetWorseAcceptParameters dwa_params;

        // How many iterations should be performed without synchronizing the current solution:
        std::uint32_t max_independent_thread_iters;

        // Number of iterations in the preliminary run - The preliminary run can be used to calibrate the parameters of the acceptance criterion
        std::uint32_t prerun_iters;

        // Maximum number of iterations
        std::uint32_t max_iters;

        // The number of consecutive iterations without improving the global best,
        // which can be considered "alarming".
        std::uint32_t iters_without_improvement_alarm;

        // The number of consecutive iterations without improving the global best,
        // after which the algorithm is terminated.
        std::uint32_t iters_without_improvement_max;

        // Maximum run time in seconds
        double max_seconds;

        // Parameters for adaptive weight adjustment
        float score_mult_accepted;
        float score_mult_improved;
        float score_mult_global_best;
        float score_decay;
        double return_to_best_known_prob;

        // Parameters that specify where to write the output
        std::string par_tuning_filename;
        std::string log_file_basename; // Will append _all.txt or _summary.txt

        Parameters() { set_defaults(); }

        explicit Parameters(const std::string& params_filename) {
            using namespace boost::property_tree;
            using ui = std::uint32_t;

            ptree t;
            read_json(params_filename, t);

            // Start out with the default parameters. Modify values using the entries in ParFileParser object
            set_defaults();

            // Generic ALNS options
            prerun_iters = t.get<ui>("prerun-iterations");
            max_iters = t.get<ui>("total-iterations");
            max_independent_thread_iters = t.get<ui>("iterations-without-syncing-threads");
            max_seconds = t.get<double>("timeout-s");
            iters_without_improvement_alarm = t.get<ui>("iters-without-improvement-alarm");
            iters_without_improvement_max = t.get<ui>("iters-without-improvement-max");

            // Acceptance criterion
            std::string ac;
            ac = t.get<std::string>("acceptance-criterion");

            if(ac == "Hill climbing") {
                acceptance_criterion_id = AcceptanceCriterionId::HillClimbing;
            } else if(ac == "Simulated annealing") {
                acceptance_criterion_id = AcceptanceCriterionId::SimulatedAnnealing;
            } else if(ac == "Threshold acceptance") {
                acceptance_criterion_id = AcceptanceCriterionId::ThresholdAcceptance;
            } else if(ac == "Great deluge") {
                acceptance_criterion_id = AcceptanceCriterionId::GreatDeluge;
            } else if(ac == "Record-to-record travel") {
                acceptance_criterion_id = AcceptanceCriterionId::RecordToRecordTravel;
            } else if(ac == "Late acceptance hill climbing") {
                acceptance_criterion_id = AcceptanceCriterionId::LateAcceptanceHillClimbing;
            } else if(ac == "Non-linear great deluge") {
                acceptance_criterion_id = AcceptanceCriterionId::NonLinearGreatDeluge;
            } else if(ac == "Worse accept") {
                acceptance_criterion_id = AcceptanceCriterionId::WorseAccept;
            } else if(ac == "Conservative worse accept") {
                acceptance_criterion_id = AcceptanceCriterionId::ConservativeWorseAccept;
            } else if(ac == "Discreet worse accept") {
                acceptance_criterion_id = AcceptanceCriterionId::DiscreetWorseAccept;
            } else if(ac == "Random walk") {
                acceptance_criterion_id = AcceptanceCriterionId::RandomWalk;
            }

            std::string pb;
            pb = t.get<std::string>("acceptance-params-base");

            if(pb == "time") {
                acceptance_params_base = AcceptanceParamsBase::Time;
            } else if(pb == "iterations") {
                acceptance_params_base = AcceptanceParamsBase::Iterations;
            }

            // Simulated annealing
            sa_params.init_accept_ratio_50p = t.get<double>("acceptance.simulated-annealing.init-ratio-50p");
            sa_params.end_accept_ratio_50p = t.get<double>("acceptance.simulated-annealing.end-ratio-50p");
            sa_params.magic_number_exponent = t.get<double>("acceptance.simulated-annealing.magic-number-exponent");
            sa_params.end_accept_ratio_refers_to_initial = t.get<bool>("acceptance.simulated-annealing.end-ratio-50p-refers-to-initial");
            sa_params.temperature_decrease_is_linear = t.get<bool>("acceptance.simulated-annealing.temperature-decrease-is-linear");
            sa_params.reheating_is_enabled = t.get<bool>("acceptance.simulated-annealing.reheating");
            sa_params.reheating_coefficient = t.get<double>("acceptance.simulated-annealing.reheating-coefficient");
            sa_params.reheating_times = t.get<double>("acceptance.simulated-annealing.reheating-times");

            // Threshold acceptance
            ta_params.start_threshold = t.get<double>("acceptance.threshold-acceptance.start-threshold");
            ta_params.end_threshold = t.get<double>("acceptance.threshold-acceptance.end-threshold");
            ta_params.threshold_decrease_is_linear = t.get<bool>("acceptance.threshold-acceptance.threshold-decrease-is-linear");

            // Great deluge
            gd_params.initial_water_level_ratio = t.get<double>("acceptance.great-deluge.initial-water-level-ratio");
            gd_params.water_level_decrease_percentage = t.get<double>("acceptance.great-deluge.water-level-decrease-pct");

            // Record to Record Travel
            rrt_params.start_deviation = t.get<double>("acceptance.record-to-record.start-deviation");
            rrt_params.end_deviation = t.get<double>("acceptance.record-to-record.end-deviation");
            rrt_params.deviation_decrease_is_linear = t.get<bool>("acceptance.record-to-record.deviation-decrease-is-linear");

            // Late acceptance hill climbing
            lahc_params.list_size = t.get<ui>("acceptance.late-acceptance-hill-climbing.list-size");
            lahc_params.accept_non_worsening = t.get<bool>("acceptance.late-acceptance-hill-climbing.allow-non-worsening");

            // Non-linear great deluge
            nlgd_params.nl_initial_water_level_ratio = t.get<double>("acceptance.non-linear-great-deluge.initial-water-level-ratio");
            nlgd_params.nl_gap_to_increase_water_level = t.get<double>("acceptance.non-linear-great-deluge.gap-to-increase-water-level");
            nlgd_params.nl_water_level_increase_pct = t.get<double>("acceptance.non-linear-great-deluge.water-level-increase-pct");
            nlgd_params.nl_water_level_decrease_exp = t.get<double>("acceptance.non-linear-great-deluge.water-level-decrease-exp-factor");

            // Worse accept
            wa_params.start_prob = t.get<double>("acceptance.worse-accept.start-probability");
            wa_params.end_prob = t.get<double>("acceptance.worse-accept.end-probability");
            wa_params.prob_decrease_is_linear = t.get<bool>("acceptance.worse-accept.probability-decrease-is-linear");

            // Conservative worse accept
            cwa_params.start_prob = t.get<double>("acceptance.conservative-worse-accept.start-probability");
            cwa_params.end_prob = t.get<double>("acceptance.conservative-worse-accept.end-probability");
            cwa_params.prob_decrease_is_linear = t.get<bool>("acceptance.conservative-worse-accept.probability-decrease-is-linear");

            // Discreet worse accept
            dwa_params.consecutive_rejects_for_100p = t.get<std::uint32_t>("acceptance.discreet-worse-accept.consecutive-rejects-for-100p");
            dwa_params.prob_increase_is_linear = t.get<bool>("acceptance.discreet-worse-accept.probability-increase-is-linear");

            // ALNS specific (heuristic selection)
            score_mult_accepted = t.get<float>("scores.accepted");
            score_mult_improved = t.get<float>("scores.improved");
            score_mult_global_best = t.get<float>("scores.global-best");
            score_decay = t.get<float>("scores.decay");

            // Back-to-best
            return_to_best_known_prob = t.get<double>("prob-return-to-best-known");

            // Output files
            par_tuning_filename = t.get<std::string>("parameter-tuning-file");
            log_file_basename = t.get<std::string>("results-log-basename");
        }

        void set_defaults() {
            sa_params =
                SimulatedAnnealingParameters(0.9,   // Initial accept ratio with 50% prob
                                             0.2,   // Final accept ratio with 50% prob
                                             1.0,   // Magic number exponent
                                             false, // Final accept ratio refers to probability vs initial solution (true)? or versus current best (false)?
                                             false, // Temperature decrease is linear (as opposed to exponential decay)
                                             false, // Enable reheating?
                                             5.0,   // Reheat to the temperature at the time of the last improvement times this coefficient
                                             3      // Reheat this number of times
                                             );

            ta_params = ThresholdParameters(0.05,  // Start Threshold
                                            0.001, // End Threshold
                                            true   // Threshold decrease is linear (as opposed to exponential decay)
                                            );

            gd_params = GreatDelugeParameters(1.1, // Initial water level ratio (* initial solution value)
                                              0.01 // Water level decrease rate (* avg gap btw current and water level)
                                              );

            rrt_params = RecordDeviationParameters(0.1,  // Start deviation
                                                   0.01, // End deviation
                                                   true  // Deviation decrease is linear (as opposed to exponential decay)
                                                   );

            lahc_params = LateAcceptanceParameters(5000, // List size
                                                   true  // Accept all non-worsening moves
                                                   );

            nlgd_params = NLGreatDelugeParameters(1.1,  // Initial water level (* initial solution value)
                                                  0.01, // Gap to increase the water level
                                                  0.05, // Water level increase percentage
                                                  0.01  // Water level decrease exponent factor
                                                  );

            wa_params = WorseAcceptParameters(0.1,  // Initial probability
                                              0.01, // Final probability
                                              true  // Probability of accepting decreases linearly (as opposed to exponentially)
                                              );

            cwa_params = ConservativeWorseAcceptParameters(0.1,  // Initial probability
                                                           0.01, // Final probability
                                                           true  // Probability of accepting decreases linearly (as opposed to exponentially)
                                                           );

            dwa_params = DiscreetWorseAcceptParameters(
                50,  // After this consecutive rejections, we accept a worsening move with probability 100%
                true // Probability of accepting increases linearly [in the number of consecutive rejects] (as opposed to exponentially)
                );

            acceptance_criterion_id = AcceptanceCriterionId::SimulatedAnnealing;
            acceptance_params_base = AcceptanceParamsBase::Iterations;
            max_iters = 50000;
            max_seconds = 60;
            iters_without_improvement_alarm = 1000;
            iters_without_improvement_max = 5000;
            score_mult_accepted = 2;
            score_mult_improved = 4;
            score_mult_global_best = 10;
            score_decay = 0.99;
            max_independent_thread_iters = 0;
            return_to_best_known_prob = -1;
            prerun_iters = 1000; // Initial iterations to "calibrate" the acceptance criterion
            par_tuning_filename = "WARNING_not_set.txt";
            log_file_basename = "WARNING_not_set";
        }
    };
}

#endif
