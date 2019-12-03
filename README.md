## Adaptive Large Neighbourhood Search

This is a header-only C++ implementation of the (Parallel) Adaptive Large Neighbourhhod Search algorithm, devised by Ropke and Pisinger.
The code is loosely based on the original implementation which Stefan Ropke has kindly shared with me.

### Usage

All parameters of the metaheuristic framework are set editing file `Params.json`.
The user needs to provide:

* An object implementing an instance of the problem, with a public method `getInstanceSize()` which returns the size of the instance (e.g., the number of vertices in a TSP instance).
* A class representing a solution of the problem, with a public method `getCost()` returning a cost to minimise (e.g., the length of a TSP tour).
* A subclass of `InitialSolutionCreator` to produce an initial solution to the problem.

The destroy and repair method will derive, respectively, from base classes `DestroyMethod` and `RepairMethod`.
The user can choose between using one the predefined acceptance criteria, or implementing his own inheriting from base class `AcceptanceCriterion`.
A handy way to monitor the algorithm and perform non-standard actions during the solution process is to provide a visitor object which inherits from `AlgorithmVisitor`.

### License

This software is licensed under the GNU Public License Version 3. See file `LICENSE` for more information.
