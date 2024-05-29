# Adaptive Large Neighbourhood Search

[![DOI](https://zenodo.org/badge/116035646.svg)](https://zenodo.org/badge/latestdoi/116035646)

This is a header-only C++ implementation of the (Parallel) Adaptive Large Neighbourhhod Search algorithm, devised by Ropke and Pisinger.
The code is loosely based on the original implementation which Stefan Ropke has kindly shared with me.

The original ALNS framework was developed in this paper:

```bib
@article{ropke_2006,
    title={An adaptive large neighborhood search heuristic for the pickup and delivery problem with time windows},
    author={Ropke, Stefan and Pisinger, David},
    year=2006,
    journal={Transportation Science},
    volume=40,
    number=4,
    pages={455--472},
    doi={10.1287/trsc.1050.0135}
}
```

This implementation also includes multiple acceptance criteria, as they were introduced in the paper [A comparison of acceptance criteria for the Adaptive Large Neighbourhood Search metaheuristic](https://link.springer.com/article/10.1007/s10732-018-9377-x) ([preprint](https://santini.in/files/papers/santini-ropke-hvattum-2017.pdf)).
If you use this software, please cite this last paper as follows:

```bib
@article{ALNS_Acceptance_Criteria,
  title={A comparison of acceptance criteria for the {Adaptive Large Neighbourhood Search} metaheuristic},
  author={Santini, Alberto and Ropke, Stefan and Hvattum, Lars Magnus},
  journal={{Journal of Heuristics}},
  volume=24,
  issue=5,
  pages={783--815},
  year=2018,
  doi={10.1007/s10732-018-9377-x}
}
```

## Usage

Functions are thoroughly documented in the code, using Doxygen syntax.
Because this is a metaheuristic framework, the user must implement the problem-specific parts.
In terms of code, this means that the user needs to provide the components described below.

### Problem instance and problem solution

The user shall implement:
* An object implementing an instance of the problem, with a public method `getInstanceSize()` which returns the size of the instance (e.g., the number of vertices in a TSP instance).
* A class representing a solution of the problem, with a public method `getCost()` returning a cost to minimise (e.g., the length of a TSP tour).
* A subclass of `InitialSolutionCreator` to produce an initial solution to the problem.

Once this is in-place, the user can create an instance of the `PALNS` object.
For example:

```cpp
auto alns = PALNS<TSPInstance, TSPSolution>{instance};
```

Where `TSPInstance` is the class modelling the problem instance and `TSPSolution` is the class modelling the problem solution.
The user must supply at least the functions mentioned above.
For example:

```cpp
std::uint32_t TSPInstance::getInstanceSize() const {
    return this->numberOfVertices;
}

double TSPSolution::getCost() const {
    return this->tourLength;
}
```

Moreover, because the ALNS method starts from a concrete solution, the user must supply an initial solution creator:

```cpp
struct TSPRandomSolutionCreator : InitialSolutionCreator<TSPSolution, TSPInstance> {
    TSPSolution create_initial_solution(const TPSInstance& instance, std::mt19937& mt) {
        return TSPSolution::shuffled_vertices(mt);
    }
}
```

For reproducibility, the initial solution creator can use a pseudo-random number generator passed by the ALNS framework, of type `std::mt19937` (from the `<random>` header).

### Destroy and repair methods

The destroy and repair methods will derive, respectively, from base classes `DestroyMethod` and `RepairMethod`.
A class inheriting from `DestroyMethod` shall implement two methods:
* `void destroy_solution(Solution sol&, std::mt19937& mt)`, which takes a solution by reference (and a PRNG, if you need it) and destroys it in-place.
* `std::unique_ptr<DestroyMethod<Solution>> clone() const`, which returns a `unique_ptr` to the base class `DestroyMethod` cloning the destroy method.

Similarly, a class inheriting from `RepairMethod` shall implement:
* `void repair_solution(Solution& sol, std::mt19937& mt)`, which takes a (destroyed) solution by reference (and a PRNG, if you need it) and repairs it in-place.
* `std::unique_ptr<RepairMethod<Solution>> clone() const`, which returns a `unique_ptr` to the base class `RepairMethod` cloning the repair method.

For example:

```cpp
struct RandomDestroy : public DestroyMethod<TSPSolution> {
    void destroy_solution(TSPSolution& solution, std::mt19937& mt) {
        std::bernoulli_distribution d(0.1);

        for(const auto& v : solution.visited_vertices()) {
            if(d(mt)) {
                solution.queue_vertex_removal(v);
            }
        }

        solution.remove_queued_vertices();
    }

    std::unique_ptr<DestroyMethod<TSPSolution>> clone() const {
        return std::make_unique<DestroyMethod<TSPSolution>>();
    }
};

struct GreedyRepair : public RepairMethod<TSPSoution> {
    void repair_solution(Solution& sol, std::mt19937& mt) {
        for(const auto& v : solution.missing_vertices()) {
            solution.insert_vertex_in_best_position(v);
        }
    }

    std::unique_ptr<RepairMethod<TSPSolution>> clone() const {
        return std::make_unique<RepairMethod<TSPSolution>>();
    }
};
```

### Other components

The user can choose between using one the predefined acceptance criteria, or implementing his own inheriting from base class `AcceptanceCriterion`.
A handy way to monitor the algorithm and perform non-standard actions during the solution process is to provide a visitor object which inherits from `AlgorithmVisitor`.
The algorithm visitor must implement the following callbacks:

* `on_algorithm_start` which is called before the first ALNS iteration.
* `on_prerun_end` which is called when the pre-run is finished (see section Parameters below).
* `on_iteration_end` which is called at the end of each iteration.
* `on_many_iters_without_improvement` which is called after a user-defined number of consecutive iterations without improving the best solution (see section Parameters below).

### Parameters

All parameters of the metaheuristic framework are set editing file `Params.json`.
They are described below:

* `prerun-iterations`: the number of initial iterations which the algorithm can use to attempt to auto-tune acceptance parameters. This is not compulsory and can be set to zero.
* `total-iterations`: a hard limit on the number of iterations.
* `timeout-s`: a hard limit on the CPU time of the algorithm.
* `iters-without-improvement-max`: a hard limit on the number of consecutive iterations without improvement to the best solution. It terminates the run early when stuck in a *flat landscape* situation.
* `iters-without-improvement-alarm`: if set, and if the user provides an algorithm visitor, after this many consecutive iterations we will call the visitor's `on_many_iters_without_improvement` method. The user can use this mechanism to try and escape local optima. For example, it can apply a stronger destroy method, or a diversification heuristic.
* `iterations-without-syncing-threads`: in the parallel version of the algorithm, this is the number of iterations during which threads do not sync information and behave like independent copies of the ALNS algorithm. Upon thread synchronisation, the accept/destroy weights and the current solution are shared.
* `prob-return-to-best-known`: at each iteration, the best known solution will replace the current solution with this probability. Zero is a good default, but the user can set a higher value to favour intensification over diversification.
* `acceptance-criterion`: one of the predefined acceptance criteria. Possible values are listed in the (unused) parameter key `acceptance-criterion-possible-values`.
* `acceptance-params-base`: some of the acceptance parameters vary during the course of the algorithm run. This parameter, which can take values `time` and `iterations` specifies if the time limit (`timeout-s`) or the iterations limit (`total-iterations`) should be considered to estimate when the run will be over.
* `scores`: these are the multipliers used to update the destroy/repair scores at each iteration. Users should provide four multipliers:
    * `global-best` if the new solution is better than the previous best solution, otherwise
    * `improved` if the new solution is better than the current solution, otherwise
    * `accepted` if the new solution was accepted by the acceptance criterion.
    * The fourth multiplier, `decay`, determines how fast the above three values affect the methods' scores. A sane default is just a little less than 1.0, to ensure that scores change smoothly.
* `parameter-tuning-file`: if performing parameter tuning, this is the filename where results will be saved.
* `results-log-basename`: basename for the result files.

#### Acceptance parameters

In the following we list parameters relative to the acceptance criteria.
We urge the users to refer to the acceptance criteria paper to better understand the role of each parameter, and the terminology used below.
The paper also present the best tuned parameters over a wide variety of problems solved.

* Simulated annealing acceptance
    * `init-ratio-50p`: solutions this much worse than the current one are accepted with 50% probability at the beginning of the run.
    * `end-ratio-50p`: solutions this much worse than the current one are accepted with 50% probability at the end of the run.
    * `end-ratio-50p-refers-to-initial`: if set to `true`, in the description of `end-ratio-50p` replace "current one" with "initial one". In other words, always evaluate the quality of new solutions compared with the initial solution and not with the current solution.
    * `temperature-decrease-is-linear`: if `true`, the decrease from the start to the end temperature is linear, otherwise it is exponential. An exponential decrease is the standard in the Simulated Annealing literature, but a linear decrease is recommended in our paper.
    * `reheating`: whether to re-increase the temperature at certain intervals.
    * `reheating-coefficient`: by how much to increase the temperature, if `reheating` is `true`.
    * `reheating-times`: how many times to increase the temperature during the run.
    * `magic-number-exponent`: the so-called "magic number" of Ropke and Pissinger's original implementation (see reference above).
* Threshold acceptance
    * `start-threshold`: threshold at the start of the run.
    * `end-threshold`: threshold at the end of the run.
    * `threshold-decrease-is-linear`: whether the threshold decreases linearly (`true`) or exponentially (`false`) between the start and end values.
* Great Deluge
    * `initial-water-level-ratio`: multiplier to the cost of the initial solution, to obtain the initial water level.
    * `water-level-decrease-pct`: percentage decrease of the water level, at each iteration.
* Record-to-record travel
    * `start-deviation`: maximum allowed deviation at the beginning of the run.
    * `end-deviation`: maximum end deviation at the end of the run.
    * `deviation-decrease-is-linear`: whether the deviation decreases linearly (`true`) or exponentially (`false`) between the start and end values.
* Late-acceptance Hill Climbing
    * `list-size`: size of the late-acceptance list.
    * `allow-non-worsening`: whether to always accept non-worsening new solution.
* Non-linear Great Deluge
    * `initial-water-level-ratio`: multiplier to the cost of the initial solution, to obtain the initial water level.
      `gap-to-increase-water-level`: threshold gap to re-increase the water level.
      `water-level-increase-pct`: percentage of the solution difference, by which the water level increases during increase phases.
      `water-level-decrease-exp-factor`: exponential factor (delta in our paper) to use when decreasing the water level.
* Worse accept
    * `start-probability`: acceptance probability at the start of the run.
    * `end-probability`: acceptance probability at the end of the run.
    * `probability-decrease-is-linear`: whether the probability decreases linearly (`true`) or exponentially (`false`) between the start and end values.

## License

This software is licensed under the GNU Public License Version 3. See file `LICENSE` for more information.
