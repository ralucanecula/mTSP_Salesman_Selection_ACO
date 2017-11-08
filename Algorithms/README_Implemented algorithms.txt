The Java source code files correspond to the implementation of six Ant Colony Optimization (ACO) algorithms investigated in our study for the
Multiple Traveling Salesman Problem (multiple-TSP).
The implementation of the investigated algorithms is based on the Java version of the ACOTSP (http://www.aco-metaheuristic.org/aco-code/) 
software package, which contains several ACO algorithms employed for solving symmetric Traveling Salesman Problem (TSP) instances.

Some of the source files have appended to their names a suffix in the form of "_nameofalgorithm", such as to uniquely identify source files
and to prevent having duplicated files with the same name.
Remove that suffix from the end of the file names, such as to be able to compile successfully the Java sources and run the algorithms.

More exactly, the meaning of "_nameofalgorithm" suffix is the following:
- "_g-MinMaxMMAS": refers to the source code for the MinMax MMAS with global-solution pheromone update with random
salesman selection (g-MinMaxMMAS) algorithm   

- "_g-MinMaxACS-d": refers to the source code for the MinMax ACS with global-solution pheromone update with deterministic
salesman selection (g-MinMaxACS-d) algorithm 

- "_g-MinMaxACS-p": refers to the source code for the MinMax ACS with global-solution pheromone update with probabilistic
salesman selection (g-MinMaxACS-p) algorithm  

- "_g-MinMaxACS-s": refers to the source code for the MinMax ACS with global-solution pheromone update with salesman and city
selected simultaneously (g-MinMaxACS-s) algorithm, and also the source code for the MinMax ACS with global-solution pheromone update with 
salesman and city selected simultaneously, with local search (g-MinMaxACS-sls) algorithm. In the latter case, the boolean ls_flag is set to
true, such that the local search procedure is applied in order to further improve the solutions constructed by ants
 
In case of source files which have no suffix in the form "_nameofalgorithm", they refer to the source code for the MinMax ACS with 
global-solution pheromone update with random salesman selection (g-MinMaxACS) algorithm.


 