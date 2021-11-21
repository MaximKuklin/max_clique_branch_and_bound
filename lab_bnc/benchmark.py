import math
import os.path
from time import time

from lab_bnc.branch_and_cut import BranchAndCutSolver
from wrapt_timeout_decorator import *

args_hard = [
    # dict(graph_path='data/DIMACS_all_ascii/C125.9.clq', iters=50),
    # dict(graph_path='data/DIMACS_all_ascii/brock200_1.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_2.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_3.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/brock200_4.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/gen200_p0.9_44.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/gen200_p0.9_55.clq', iters=50),
    dict(graph_path='data/DIMACS_all_ascii/p_hat1000-1.clq', iters=50),
    # dict(graph_path='data/DIMACS_all_ascii/san1000.clq', iters=50)
]



def main():

    @timeout(5000)
    def get_solution():
        start = time()
        bnb.cut()
        end = time()
        return start, end

    if os.path.exists("results.txt"):
        raise FileExistsError

    file = open("results.txt", 'w')
    file.write("CPU = AMD Ryzen 7 2700X Eight-Core Processor, CPU MHz = 2200, Time Limit = 3600 sec.\n")
    file.write("Graph Time Clique_Size\n")

    for argument in args_hard:

        name = os.path.basename(argument['graph_path'])

        full_start = time()

        bnb = BranchAndCutSolver(mode="LP", **argument)

        try:
            start, end = get_solution()
        except:
            print("Oops, timeout!!!!!!!!!!!!!!!!!!!")
            continue

        print(f"Final max clique has size {math.floor(bnb.best_solution)}")
        print(f"Vertexes are: {bnb.best_vertexes}")

        print(f"\nFull Time {end - full_start:.4f} sec.")
        print(f"CPLEX Time {end - start:.4f} sec.")

        final_string = ' '.join([name, f'{end - full_start:.3f}', f"{math.floor(bnb.best_solution)}"])

        file.write(final_string + "\n")

    file.close()


if __name__ == '__main__':
    main()