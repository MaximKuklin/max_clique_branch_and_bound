# Solve Max Clique Problem using CPLEX and Branch & Bound Method


## Run code

To run BnB method on some DIMASC graph, use the following command:

`python3 optimization/branch_and_bound.py -p data/DIMACS_all_ascii/MANN_a9.clq`

## Logs

After running the command above, you will see the output results like this:

```text
Best founded clique using heuristic has size 16
Final max clique has size 16
Vertexes are: [33, 36, 4, 5, 39, 7, 9, 42, 8, 12, 15, 18, 21, 24, 27, 30]

Full Time 0.2593 sec.
CPLEX Time 0.1016 sec.
```


## Result table

|       Graph       | BnB time (sec) | Full time (sec) | Heurictic clique size | BnB max clique size | True max clique size | Coloring iters |
|:-----------------:|:--------------:|:---------------:|:---------------------:|:-------------------:|:--------------------:|:--------------:|
|  johnson8-2-4.clq |     0.0044     |       0.17      |           4           |          4          |           4          |        5       |
| johnson16-2-4.clq |     0.0037     |       1.38      |           8           |          8          |           8          |       16       |
|   hamming8-4.clq  |      18.18     |      35.48      |           16          |          16         |          16          |       30       |
|   c-fat200-1.clq  |     0.0455     |       4.09      |           12          |          12         |           -          |        1       |
|   c-fat200-2.clq  |     0.0469     |       2.01      |           24          |          24         |           -          |        1       |
|   c-fat200-5.clq  |      4.20      |      10.32      |           58          |          58         |           -          |       10       |
|   c-fat500-1.clq  |      0.302     |      116.22     |           14          |          14         |           -          |        1       |
|   c-fat500-2.clq  |      0.320     |       60.3      |           26          |          26         |           -          |        1       |
|    MANN_a9.clq    |     0.1056     |       0.26      |           16          |          16         |          16          |        1       |
|    keller4.clq    |     118.44     |      128.44     |           8           |          11         |          11          |       50       |

