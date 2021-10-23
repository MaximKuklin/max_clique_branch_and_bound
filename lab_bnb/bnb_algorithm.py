import networkx as nx
from typing import Union
import cplex
from copy import deepcopy
import itertools
from igraph import Graph
import numpy as np
from time import time

from optimization.lab_bnb.lab1 import MaxCliqueSolver

EPS = 1e-5


class BranchAndBoundSolver(MaxCliqueSolver):
    def __init__(self, mode, graph_path, sense='maximize'):
        super().__init__(mode, graph_path, sense)

        self.solution = None
        self.best_solution = 0
        self.best_vertexes = None
        self.branch_num = 0

        self._first_step()

    def update_constraints(self):
        pass

    def _first_step(self):
        self.clique_problem = MaxCliqueSolver(mode=self.mode, graph_path=self.path)
        names, obj, lower_bounds, upper_bounds = self.clique_problem.get_variables()
        self.set_variables(names, obj, lower_bounds, upper_bounds)

        constraints = self.get_constraints()
        self.set_constraints(constraints, "L", rhos=1)

        self.create_problem()
        # solution, t = self.solve()

        # variables = solution.get_values()

        # self.variables = variables

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
            max_idx = self.variables.index(max(float_vars))
            return max_idx
        else:
            return None

    def branching_largest_first(self):

        solution = self.solve()[0]
        self.variables = solution.get_values()

        if int(sum(self.variables) + EPS) > self.best_solution:

            max_idx = self._get_branch_variable()

            if max_idx is None:
                self.best_solution = sum(self.variables)
                self.best_vertexes = self.variables
                return self.best_solution
            else:
                self.branch_num += 1
                current_branch = self.branch_num

                self._create_branch_constraint(max_idx, 1.0, current_branch)
                branch_1 = self.branching_largest_first()
                self.problem.linear_constraints.delete(f'branch_{current_branch}')

                self._create_branch_constraint(max_idx, 0.0, current_branch)
                branch_2 = self.branching_largest_first()
                self.problem.linear_constraints.delete(f'branch_{current_branch}')
                return max(branch_1, branch_2)
        return 0


bnb = BranchAndBoundSolver(mode="LP", graph_path="data/DIMACS_all_ascii/debug.clq")

bnb.branching_largest_first()

# print(answer)
print(bnb.best_solution)
print(bnb.best_vertexes)
