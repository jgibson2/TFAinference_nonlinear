import csv
import time

from mystic.solvers import LatticeSolver, PowellDirectionalSolver

from mysticRealMain import *

"""ProcessPool is a way that the optimization jobs can be split up for parallelization, ParallelPool is an alternative to ProcessPool, one or the other may be faster depending on the job"""
from pathos.multiprocessing import ProcessPool as Pool

# from pathos.pools import ParallelPool as Pool

"""imports matrices as .csv files and makes lists of lists of numbers to represent them
binaryCS matrix dimensions are numGenes x numTFs
binaryTFA matrix dimensions are numTFs x numSamples
geneExpression matrix dimensions are numGenes x numSamples
All .csv files are assumed to only have numbers and no labels
If .csv files have labels, change 'for i in range(len(your_list))' to 'for i in range(1,len(your_list))' and change 'for j in range(len(your_list[0]))' to 'for j in range(1,len(your_list[0]))'"""
with open("Data//binaryCS.csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
signedBinaryCS = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(int(your_list[i][j]))
    signedBinaryCS.append(row)

with open("Data//binaryTFAtrimmed.csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
binaryTFA = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(int(your_list[i][j]))
    binaryTFA.append(row)

with open("Data//geneExpressiontrimmed.csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
geneExpressionMatrix = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        try:
            row.append(float(your_list[i][j]))
        except:
            for k in range(len(your_list[i][j])):
                if your_list[i][j][k] == "e" or your_list[i][j][k] == "E":
                    row.append(float(your_list[i][j][:k]) * (10 ** float(your_list[i][j][k + 1:])))
    geneExpressionMatrix.append(row)

numGenes = len(signedBinaryCS)
numTFs = len(signedBinaryCS[0])
numSamples = len(binaryTFA[0])

"""createParams is used fro the original non-linear model, if the updated non-linear model is used, change createParams to createParamsUpdate1, if the bi-linear model is used, change createParams to createParamsBL"""
params = createParams(signedBinaryCS, binaryTFA)
dim = len(params)

proc = process(signedBinaryCS, binaryTFA, geneExpressionMatrix, params)

"""sets the upper and lower bounds to be lists of length equal to the number of parameters"""
(lb, ub) = createBounds(params)

my_solver = PowellDirectionalSolver(dim)  # initializes the PowellDirectionalSolver with the number of parameters
my_solver.SetStrictRanges(lb, ub)  # sets the range that the solutions must be within for PowellDirectionalSolver
start = time.time()
solver = LatticeSolver(dim,
                       nbins=2)  # initalizes the LatticeSolver with the number of parameters and the number of bins for starting points
solver.SetStrictRanges(lb, ub)  # sets the range that the solutions must be within for the LatticeSolver
solver.SetNestedSolver(my_solver)  # nests the PowellDirectionalSolver within the LatticeSolver
"""sets a mapper which breaks the optimization into various distinct tasks which can be run independently 
Sets the number of nodes available to the solver equal to the number of nodes * the number of cpus per node"""
solver.SetMapper(Pool(nodes=200).map)
"""processDataForOptimizing is used for the original non-linear model, if the updated non-linear model is used, change processDataForOptimizing to processDataForOptimizingUpdate1, if the bi-linear model is used, change processDataForOptimizing to processDataForOptimizingBL"""
solver.Solve(proc.processDataForOptimizing)  # LatticeSolver optimizes on the function processDataForOptimizing
sol = solver.Solution()  # sets sol to the list of values that each parameter is optimized to
stop = time.time()
lst = []
for j in range(len(sol)):
    lst.append(params[j] + "->" + str(
        sol[j]))  # the solution is just a list of values, so this assigns each parameter to its corresponding value
lst.append("time->" + str((stop - start) / 3600))  # the time in hours that the optimization took
lst.append("error->" + str(proc.processDataForOptimizing(sol)))  # the error that the optimizer got
with open("paramRulesMystic//paramRulesReal",
          "w") as file:  # saves the parameters with their values, time, and error as a file
    file.write(str(lst))
