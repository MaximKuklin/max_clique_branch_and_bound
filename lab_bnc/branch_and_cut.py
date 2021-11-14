import argparse
import itertools
import math

import networkx as nx
import numpy as np
from time import time
from lab_bnb.max_clique_solver import MaxCliqueSolver
from lab_bnb.branch_and_bound import parse_args, BranchAndBoundSolver

EPS = 1e-5

import random
random.seed(0)        # or any integer
np.random.seed(0)


class BranchAndCutSolver(MaxCliqueSolver):
    def __init__(self, mode, graph_path, iters, sense='maximize'):
        super().__init__(mode, graph_path, iters, sense)

        self.solution = None
        self.best_solution = 0
        self.best_vertexes = None
        self.branch_num = 0
        self.cut_num = 0
        self.total_time = 0
        self.max_time = 3600  # seconds

        self._first_step()

    def _first_step(self):

        names, obj, lower_bounds, upper_bounds = super().get_variables()
        self.set_variables(names, obj, lower_bounds, upper_bounds)

        constraints = self.get_constraints()
        self.set_constraints(constraints, "L", rhos=1)

        self.create_problem()
        self.set_complement()
        self.complement_graph = nx.Graph(self.complement)

    def get_constraints(self):

        constraints = []

        independent_constraints = \
            [self._get_constraint(ind_set, self.names, use_ind=False) for ind_set in self.ind_sets]
        single_constraints = \
            [self._get_constraint([single_vertex], self.names, use_ind=False) for single_vertex in range(len(self.names))]

        constraints.extend(independent_constraints)
        constraints.extend(single_constraints)
        return independent_constraints

    def _create_branch_constraint(self, idx, constr, curr_branch_num):

        constraint = [[idx], [1.0]]

        self.problem.linear_constraints.add(
            lin_expr=[constraint],
            senses=['E'],
            rhs=[constr],
            names=[f"branch_{curr_branch_num}"]
        )

    def add_constraint(self, vertexes: list):
        n = len(vertexes)
        names = [f"x{i}" for i in vertexes]
        constraint = [names, [1.0]*n]

        self.problem.linear_constraints.add(
            lin_expr=[constraint],
            senses=['L'],
            rhs=[1.0],
            names=[f"cut_{self.cut_num}"]
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

    def cut(self):
        solution, t = self.solve()

        # self.total_time += t
        # if self.total_time >= self.max_time:
        #     print("!!! Time limit reached, turning back !!!")
        #     return

        self.variables = solution.get_values()

        if int(sum(self.variables) + EPS) <= self.best_solution:
            return

        while True:

            result = self.separation(self.variables)

            if result is None:
                break
            else:
                self.add_constraint(result)
                self.cut_num += 1

            solution, t = self.solve()
            self.variables = solution.get_values()

            if int(sum(self.variables) + EPS) <= self.best_solution:
                return

        self.branching()

    # def check_solution(self, nodes):
    #     def separation(i, j):
    #         independent_vertexes = set([i for i in range(len(self.variables))])
    #         new_constraints = [i, j]
    #         for vertex in (i, j):
    #             non_neighbors = set(self.complement_graph.neighbors(vertex))
    #             independent_vertexes.intersection_update(non_neighbors)
    #
    #     for i, j in itertools.combinations(nodes, 2):
    #         if not self.graph.has_edge(i, j):
    #             self.add_constraint()
    #             return False
    #     return True

    def branching(self):
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
            self.branching()
            self.problem.linear_constraints.delete(f'branch_{current_branch}')

            self._create_branch_constraint(idx, 0.0, current_branch)
            self.branching()
            self.problem.linear_constraints.delete(f'branch_{current_branch}')

    def separation(self, values):
        # heuristic to find independent set with max weight
        values = np.array(values)

        current_vertex = np.argmax(values)
        weights_sum = values[current_vertex]
        independent_vertexes = set([i for i in range(len(values))])
        new_constraints = [current_vertex]
        while True:
            non_neighbors = set(self.complement_graph.neighbors(current_vertex))
            independent_vertexes.intersection_update(non_neighbors)

            weights = values[list(independent_vertexes)]
            if np.any(weights > 0):
                max_idx = np.argmax(weights)
                current_vertex = list(independent_vertexes)[max_idx]
                new_constraints.append(current_vertex)
                weights_sum += weights[max_idx]
            else:
                if len(new_constraints) > 2 and weights_sum > 1:
                    return new_constraints
                else:
                    return None


def main():
    args = parse_args()
    bnc = BranchAndCutSolver(mode="LP", graph_path=args.path, iters=args.iter_coloring)
    bnc.cut()


if __name__ == '__main__':
    main()