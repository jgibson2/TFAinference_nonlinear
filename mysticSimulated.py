import csv
import sys
import time

from mystic.solvers import LatticeSolver, PowellDirectionalSolver

from mysticSimulatedMain import *

"""ProcessPool is a way that the optimization jobs can be split up for parallelizaton, ParallelPool is an alternative to ProcessPool, one or the other may be faster depending on the job"""
from pathos.multiprocessing import ProcessPool as Pool

# from pathos.pools import ParallelPool as Pool


"""sets ide to the sbatch array number 
ide is only used for importing different simulated files and creating different simulated data"""
ide = int(sys.argv[1])

"""imports matrices and vectors as .csv files and makes lists of lists of numbers to represent matrices and makes lists of numbers to represent vectors
binaryCS matrix dimensions are numGenes x numTFs
binaryTFA matrix dimensions are numTFs x numSamples
simulatedControlStrength matrix dimensions are numGenes x numTFs
simulatedActivity matrix dimensions are numTFs x numSamples
simulatedBaselines vector is length numGenes
simulatedScaling vector is length numGenes
simulatedNoise matrix dimensions are numGenes x numSamples
All .csv files are assumed to only have numbers and no lables
If .csv files have labels, change 'for i in range(len(your_list))' to 'for i in range(1,len(your_list))' and change 'for j in range(len(your_list[0]))' to for j in range(1,len(your_list[0]))'"""
with open("SimulatedData//binaryCS.csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
binaryCS = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(int(your_list[i][j]))
    binaryCS.append(row)

"""only used for importing and creating different simulated data
Determines which binaryTFA gets imported"""
if (ide - 1) % 9 < 3:
    numSamp = 12
elif (ide - 1) % 9 > 2 and (ide - 1) % 9 < 6:
    numSamp = 24
else:
    numSamp = 36

with open("SimulatedData//binaryTFA" + str(numSamp) + "Samples.csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
binaryTFA = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(int(your_list[i][j]))
    binaryTFA.append(row)

numGenes = len(binaryCS)
numTFs = len(binaryCS[0])
numSamples = len(binaryTFA[0])

with open("SimulatedData//SimulatedControlStrength" + str((ide - 1) % 54 + 1) + ".csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
simulatedControlStrength = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(float(your_list[i][j]))
    simulatedControlStrength.append(row)
signedBinaryCS = sign(simulatedControlStrength)

with open("SimulatedData//SimulatedActivity" + str((ide - 1) % 54 + 1) + ".csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
simulatedActivity = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(float(your_list[i][j]))
    simulatedActivity.append(row)

with open("SimulatedData//SimulatedBaselines" + str((ide - 1) % 54 + 1) + ".csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
simulatedBaselines = []
for i in range(len(your_list)):
    simulatedBaselines.append(float(your_list[i][0]))

with open("SimulatedData//SimulatedScaling" + str((ide - 1) % 54 + 1) + ".csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
simulatedScaling = []
for i in range(len(your_list)):
    simulatedScaling.append(float(your_list[i][0]))

with open("SimulatedData//SimulatedNoise" + str((ide - 1) % 54 + 1) + ".csv") as f:
    reader = csv.reader(f)
    your_list = list(reader)
simulatedNoise = []
for i in range(len(your_list)):
    row = []
    for j in range(len(your_list[0])):
        row.append(float(your_list[i][j]))
    simulatedNoise.append(row)

"""calculates the simulated expression matrix with either the original non-linear model or the updated non-linear model
Only used for generating different simulated data"""
if (ide - 1) % 18 < 9:
    simulatedExpressionMatrix = calculateExpression(simulatedControlStrength, simulatedActivity, simulatedBaselines,
                                                    simulatedScaling)
else:
    simulatedExpressionMatrix = calculatedExpression1(simulatedControlStrength, simulatedActivity, simulatedBaselines,
                                                      simulatedScaling)
for i in range(len(simulatedExpressionMatrix)):
    for j in range(len(simulatedExpressionMatrix[1])):
        simulatedExpressionMatrix[i][j] += simulatedNoise[i][j]

"""the error expression is calculated with either the original non-linear model, the updated non-linear model, or the bi-linear model,
so cP is assigneds one of the 3 functions so that it will create the list of parameters with the correct model"""
if ide < 55:
    cP = createParams
elif ide > 54 and ide < 109:
    cP = createParamsUpdate1
else:
    cP = createParamsBL

params = cP(signedBinaryCS, binaryTFA)
dim = len(params)
(lb, ub) = createBounds(params)

proc = process(signedBinaryCS, binaryTFA, simulatedExpressionMatrix, params)

"""assigns g to be one of the 3 processDataForOptimizing functions so that Mystic optimizes on the correct one
Only used for generating different simulated data"""
if ide < 55:
    g = proc.processDataForOptimizing
elif ide > 54 and ide < 109:
    g = proc.processDataForOptimizingUpdate1
else:
    g = proc.processDataForOptimizingBL

my_solver = PowellDirectionalSolver(dim)  # initializes the PowellDirectionalSolver with the number of parameters
my_solver.SetStrictRanges(lb, ub)  # sets the range that the solutions must be within for PowellDirectionalSolver
start = time.time()
solver = LatticeSolver(dim,
                       nbins=2)  # initalizes the LatticeSolver with the number of parameters and the number of bins for starting points
solver.SetStrictRanges(lb, ub)  # sets the range that the solutions must be within for the LatticeSolver
solver.SetNestedSolver(my_solver)  # nests the PowellDirectionalSolver within the LatticeSolver
solver.SetMapper(Pool(nodes=100).map)
solver.Solve(g)  # LatticeSolver optimizes on the function g
sol = solver.Solution()  # sets sol to the list of values that each parameter is optimized to
stop = time.time()
lst = []
for j in range(len(sol)):
    lst.append(params[j] + "->" + str(
        sol[j]))  # the solution is just a list of values, so this assigns each parameter to its corresponding value
lst.append("time->" + str((stop - start) / 3600))  # the time in hours that the optimization took
lst.append("error->" + str(g(sol)))  # the error that the optimizer got
with open("paramRulesMystic//paramRulesSimulated" + str(ide),
          "w") as file:  # saves the parameters with their values, time, and error as a file
    file.write(str(lst))
