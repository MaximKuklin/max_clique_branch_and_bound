import math
import os.path
from time import time

from optimization.lab_bnb.branch_and_bound import BranchAndBoundSolver


args = [
    dict(graph_path='data/DIMACS_all_ascii/johnson8-2-4.clq', iters=5),
    dict(graph_path='data/DIMACS_all_ascii/johnson16-2-4.clq', iters=25),
    dict(graph_path='data/DIMACS_all_ascii/hamming8-4.clq', iters=30),
    dict(graph_path='data/DIMACS_all_ascii/c-fat200-1.clq', iters=1),
    dict(graph_path='data/DIMACS_all_ascii/c-fat200-2.clq', iters=1),
    dict(graph_path='data/DIMACS_all_ascii/c-fat200-5.clq', iters=10),
    dict(graph_path='data/DIMACS_all_ascii/c-fat500-1.clq', iters=1),
    dict(graph_path='data/DIMACS_all_ascii/c-fat500-2.clq', iters=1),
    dict(graph_path='data/DIMACS_all_ascii/MANN_a9.clq', iters=1),
    dict(graph_path='data/DIMACS_all_ascii/keller4.clq', iters=50),
]

args_hard = [
    dict(graph_path='data/DIMACS_all_ascii/C125.9.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_1.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_2.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_3.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_4.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/gen200_p0.9_44.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/gen200_p0.9_55.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/p_hat1000-1.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/san1000.clq', iters=50)
]

def main():

    if os.path.exists("results.txt"):
        raise FileExistsError

    file = open("results.txt", 'w')
    file.write("CPU = AMD Ryzen 7 2700X Eight-Core Processor, CPU MHz = 2200, Time Limit = 3600 sec.\n")
    file.write("Graph Time Clique_Size\n")

    for argument in args:

        name = os.path.basename(argument['graph_path'])

        full_start = time()
        bnb = BranchAndBoundSolver(mode="LP", **argument)

        start = time()
        bnb.branching_largest_first()
        end = time()

        print(f"Final max clique has size {math.floor(bnb.best_solution)}")
        print(f"Vertexes are: {bnb.best_vertexes}")

        print(f"\nFull Time {end - full_start:.4f} sec.")
        print(f"CPLEX Time {end - start:.4f} sec.")

        final_string = ' '.join([name, f'{end - full_start:.3f}', f"{math.floor(bnb.best_solution)}"])

        file.write(final_string + "\n")

    file.close()


if __name__ == '__main__':
    main()