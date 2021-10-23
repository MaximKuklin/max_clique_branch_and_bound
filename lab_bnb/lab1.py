import argparse
import os
import networkx as nx
from typing import Union
import cplex
from copy import deepcopy
import itertools
from igraph import Graph
import numpy as np
from time import time

EPS = 1e-5


class Solver:
    def __init__(self, mode, sense='maximize'):

        assert mode in ["ILP", "LP"]
        self.mode = mode
        self.sense = sense
        self.problem = None

        self.names = []
        self.obj = []
        self.lower_bounds = []
        self.upper_bounds = []

    def create_problem(self):
        self.problem = cplex.Cplex()

        if self.sense == 'maximize':
            self.problem.objective.set_sense(self.problem.objective.sense.maximize)
        elif self.sense == 'minimize':
            self.problem.objective.set_sense(self.problem.objective.sense.minimize)

        self.problem.variables.add(
            obj=self.obj,
            lb=self.lower_bounds,
            ub=self.upper_bounds,
            names=self.names
        )

        self.problem.linear_constraints.add(
            lin_expr=self.constraints,
            senses=self.constraint_senses,
            rhs=self.rhos,
            names=self.constraint_names
        )

        self._set_mode()

        self.problem.set_log_stream(None)
        self.problem.set_results_stream(None)

    def set_variables(self, names: list, obj_values: list, lower_bounds: list, upper_bounds: list):
        self.names = names
        self.obj = obj_values
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds

    def set_constraints(self,
                        constraints: list, constraint_senses: Union[list, str],
                        rhos: Union[list, int], constraint_names: Union[list] = None):

        self.constraints = constraints

        if isinstance(constraint_senses, str):
            self.constraint_senses = [constraint_senses] * len(constraints)

        if isinstance(rhos, int):
            self.rhos = [rhos] * len(constraints)

        if constraint_names is None:
            self.constraint_names = [f"с{i}" for i in range(0, len(self.constraints))]

    def _set_mode(self):
        if self.mode == "ILP":
            for i in range(len(self.obj)):
                self.problem.variables.set_types(i, self.problem.variables.type.binary)
        if self.mode == "LP":
            for i in range(len(self.obj)):
                self.problem.variables.set_types(i, self.problem.variables.type.continuous)

        for i in range(len(self.rhos)):
            if self.mode == "ILP":
                self.rhos[i] = int(self.rhos[i])
            elif self.mode == "LP":
                self.rhos[i] = float(self.rhos[i])

        assert len(self.obj) == len(self.upper_bounds) == len(self.lower_bounds)

        if self.mode == "ILP":
            for i in range(len(self.obj)):
                self.obj[i], self.upper_bounds[i], self.lower_bounds[i] = \
                    int(self.obj[i]), int(self.upper_bounds[i]), int(self.lower_bounds[i])
        elif self.mode == "LP":
            for i in range(len(self.obj)):
                self.obj[i], self.upper_bounds[i], self.lower_bounds[i] = \
                    float(self.obj[i]), float(self.upper_bounds[i]), float(self.lower_bounds[i])

    def solve(self):
        start = time()
        self.problem.solve()
        end = time()
        solution = self.problem.solution

        t = end - start
        return solution, t


class MaxCliqueSolver(Solver):
    def __init__(self, mode, graph_path, sense='maximize'):
        super().__init__(mode, sense)

        self.path = graph_path
        self.n_vertex = 0
        self.n_edges = 0

        self.graph = self.get_graph()
        self.ind_sets = []
        self.get_ind_sets()

    def get_graph(self):
        assert self.path.endswith(".clq")

        with open(self.path, 'r') as graph:
            data = graph.readlines()

        edges = []
        for line in data:
            if line.startswith("p"):
                graph_info = line[:-1].split(' ')
                n, m = graph_info[-2:]
                n, m = int(n), int(m)
            elif line.startswith("e"):
                edge = line[:-1].split(' ')[-2:]
                edge = (int(edge[0]) - 1, int(edge[1]) - 1)  # to start from 0
                edges.append(edge)
            elif line.startswith('c'):
                print(line)

        g = nx.Graph(edges)

        self.graph = g
        self.n_vertex = n
        self.n_edges = m

        return g

    def get_ind_sets(self):
        strategies = [nx.coloring.strategy_largest_first,
                      nx.coloring.strategy_random_sequential,
                      nx.coloring.strategy_independent_set,
                      nx.coloring.strategy_connected_sequential_bfs,
                      nx.coloring.strategy_connected_sequential_dfs,
                      nx.coloring.strategy_saturation_largest_first]

        for strategy in strategies:
            d = nx.coloring.greedy_color(self.graph, strategy=strategy)
            for color in set(color for node, color in d.items()):
                self.ind_sets.append(
                    [key for key, value in d.items() if value == color])

    def _get_constraint(self, ind, names, use_ind=True):
        size = len(names)
        if use_ind:
            indexes = list(range(size))
        else:
            indexes = deepcopy(names)
        values = [0] * size
        for i in ind:
            values[i] = 1
        return [indexes, values]

    def get_variables(self):
        names = [f"x{i}" for i in range(0, self.n_vertex)]
        obj, lower_bounds, upper_bounds = [1]*self.n_vertex, [0]*self.n_vertex, [1]*self.n_vertex
        return names, obj, lower_bounds, upper_bounds

    def get_constraints(self):
        inverted = list(nx.complement(self.graph).edges())

        f_e = inverted[0]
        first_constraint = self._get_constraint(f_e, self.names, use_ind=False)
        other_constraints = [self._get_constraint(edge, self.names) for edge in inverted[1:]]
        independent_constraints = [self._get_constraint(ind_set, self.names) for ind_set in self.ind_sets]

        constraints = [first_constraint]
        constraints.extend(other_constraints)
        constraints.extend(independent_constraints)

        return constraints

def main():
    start_full = time()
    clique_problem = MaxCliqueSolver(mode="ILP", graph_path="data/DIMACS_all_ascii/brock200_2.clq")

    names, obj, lower_bounds, upper_bounds = clique_problem.get_variables()
    clique_problem.set_variables(names, obj, lower_bounds, upper_bounds)

    constraints = clique_problem.get_constraints()
    clique_problem.set_constraints(constraints, "L", rhos=1)

    clique_problem.create_problem()
    solution, t = clique_problem.solve()
    end_full = time()

    clique = np.argwhere(np.array(solution.get_values()) > 1 - EPS).squeeze()
    print(f"Clique vertexes:\n{clique}\n")
    print(f"Max clique size: {clique.size}\n")

    print(f"{end_full - start_full:.4f} s. was spent on full pipeline")
    print(f"{t:.4f} s. was spent on solver")

if __name__ == '__main__':
    main()