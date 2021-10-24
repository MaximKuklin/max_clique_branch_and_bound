import argparse
import math
import numpy as np
from time import time

from optimization.lab_bnb.max_clique_solver import MaxCliqueSolver

EPS = 1e-5


class BranchAndBoundSolver(MaxCliqueSolver):
    def __init__(self, mode, graph_path, iters, sense='maximize'):
        super().__init__(mode, graph_path, iters, sense)

        self.solution = None
        self.best_solution = 0
        self.best_vertexes = None
        self.branch_num = 0
        self.total_time = 0
        self.max_time = 3600  # seconds

        self._first_step()

    def update_constraints(self):
        pass

    def _first_step(self):

        names, obj, lower_bounds, upper_bounds = super().get_variables()
        self.set_variables(names, obj, lower_bounds, upper_bounds)

        constraints = self.get_constraints()
        self.set_constraints(constraints, "L", rhos=1)

        self.create_problem()

        self.best_vertexes, self.best_solution = self.greedy_clique_heuristic()

        print(f"Best founded clique using heuristic has size {self.best_solution}")


    def _create_branch_constraint(self, idx, constr, curr_branch_num):

        constraint = [[idx], [1.0]]

        self.problem.linear_constraints.add(
            lin_expr=[constraint],
            senses=['E'],
            rhs=[constr],
            names=[f"branch_{curr_branch_num}"]
        )

    def _get_branch_variable(self):
        float_vars = list(filter(lambda x: 0+EPS < x < 1-EPS, self.variables))

        if float_vars:
            closest_to_zero_fn = (lambda list_value: abs(list_value - 0))
            closest_to_one_fn = (lambda list_value: abs(list_value - 1))

            closest_to_zero = min(float_vars, key=closest_to_zero_fn)
            closest_to_one = 1 - min(float_vars, key=closest_to_one_fn)

            if closest_to_zero < 1 - closest_to_one:
                closest = closest_to_zero
            else:
                closest = closest_to_one

            idx = self.variables.index(closest)
            return idx
        else:
            return None

    def branching_largest_first(self):

        solution, t = self.solve()
        self.total_time += t
        if self.total_time >= self.max_time:
            print("!!! Time limit reached, turning back !!!")
            return

        self.variables = solution.get_values()

        if int(sum(self.variables) + EPS) <= self.best_solution:
            return

        idx = self._get_branch_variable()

        if idx is None:
            curr_clique = np.argwhere(
                (1-EPS < np.array(self.variables)) & (np.array(self.variables) < 1+EPS)
            ).squeeze().tolist()

            if self.is_clique(curr_clique):
                self.best_solution = sum(self.variables)
                self.best_vertexes = curr_clique
                print(f"Current best solution: {math.floor(self.best_solution)}")
        else:
            self.branch_num += 1
            current_branch = self.branch_num

            self._create_branch_constraint(idx, 1.0, current_branch)
            self.branching_largest_first()
            self.problem.linear_constraints.delete(f'branch_{current_branch}')

            self._create_branch_constraint(idx, 0.0, current_branch)
            self.branching_largest_first()
            self.problem.linear_constraints.delete(f'branch_{current_branch}')

        return


def parse_args():
    args = argparse.ArgumentParser('Solve Max Clique Problem using BnB')
    args.add_argument('-p', '--path', required=True, help='Path to file DIMACS file')
    args.add_argument('--iter_coloring', required=False, default=50, type=int,
                      help="Set how many times run graph coloring algorithm")
    args = args.parse_args()

    return args


def main():

    args = parse_args()

    full_start = time()
    bnb = BranchAndBoundSolver(mode="LP", graph_path=args.path, iters=args.iter_coloring)

    start = time()
    bnb.branching_largest_first()
    end = time()

    print(f"Final max clique has size {math.floor(bnb.best_solution)}")
    print(f"Vertexes are: {bnb.best_vertexes}")

    print(f"\nFull Time {end - full_start:.4f} sec.")
    print(f"CPLEX Time {end - start:.4f} sec.")


if __name__ == '__main__':
    main()
