Andrew Krause
8/1/2018
akrause@mit.edu

Prerequisites to run: need local installation of mystic, numpy, time, csv, pathos, python 2.7.12

To run on real data:

mysticReal.py -- takes in gene expression matrix, signed binary CS matrix, and binary TFA matrix and runs mystic optimization on them to generate param rules
saves the list of param rules to a file named paramRulesReal in a folder named paramRulesMystic which is in the same directory as mysticReal.py
imports a class and functions from mysticRealMain.py
gene expression matrix, signed binary CS matrix, and binary TFA matrix should be saved as .CSV files in a directory that is in same directory as mysticReal.py
all .CSV files are imported into python and represented as lists of numbers or lists of lists of numbers
binaryCS dimensions should be numGenes x numTFs, binaryTFA dimensions should be numTFs x numSamples, geneExpressionMatrix dimensions should be numGenes x numSamples 
number of nodes set in SetMapper can be changed, more nodes will run faster until a max number of nodes required is reach
number of nodes given to SetMapper should equal number of nodes * number of cpus per node which are set in mysticReal.sbatch
increasing number of bins given to LatticeSolver will increase accuracy but will also increase time
running a network of 278 genes, 34 TFs, and 44 samples with 4 bins took about 14 hours

mysticRealMain.py -- defines functions transpose, createParams, createParamsUpdate1, createParamsBL, createBounds, and class process with functions processDataForOptimizing, processDataForOptimizingUpdate1, and processDataForOptimizingBL
tranpose takes a list of lists of numbers as input and outputs a tranposed list of lists of numbers
createParams, createParamsUpdate1, and createParamsBL take signedBinaryCS and binaryTFA as input and output of list of strings that represent the parameters that are optimized
createBounds takes params (a list of strings) as input and outputs two lists, one for the lower bound and one for the upper bounds
process takes signedBinaryCS, binaryTFA, geneExpressionMatrix, and params as input
processDataForOptimizing, processDataForOptimizingUpdate1, and processDataForOptimizingBL all take a list of numbers of length equal to the number of params as input

mysticReal.sbatch -- runs mysticReal.py, to execute on cluster "sbatch mysticReal.sbatch"
value assigned to nodes * value assigned to cpus per node should equal the number of nodes set in SetMapper in mysticReal.py

fixParamRulesReal.py -- takes in paramRulesReal output and converts it into a list of param rules that Mathematica can use
saves the fixed list of rules to a file named fixParamRulesReal in the folder paramRulesMystic

fixParamRulesReal.sbatch -- runs fixParamRulesReal.py, to execute on cluster "sbatch fixParamRulesReal.sbatch"



To generate and run on simulated data:

mysticSimulated.py -- takes in binaryCS, binaryTFA, simulatedControlStrength, simulatedActivity, and simulatedNoise matrices and simulatedBaselines and simulatedScaling vectors and runs mystic optimization on them to generate param rules
saves the list of param rules to a file named mysticSimulated in a folder named paramRulesMystic which is in the same directory as mysticSimulated.py
imports a class and functions from mysticRealMain.py
binaryCS, binaryTFA, simulatedControlStrength, simulatedActivity, simulatedNoise, simulatedBaselines, and simulatedScaling should be saved as .CSV files in a director that is in the same directory as mysticReal.py
all .CSV files are imported into python and represented as lists of numbers or lists of lists of numbers
binaryCS dimensions are numGenes x numTFs, binaryTFA dimensions are numTFs x numSamples, simulatedControlStrength dimensions are numGenes x numTFs, simulatedActivity dimensions are numTFs x numSamples, simulatedNoise dimensions are numGenes x numSamples, simulatedBaselines has length numGenes, simulatedScaling has length numGenes
ide is used to import and generate different sets of simulated data
if (ide - 1) % 18 < 9, calculateExpression is used to generate the simulatedExpressionMatrix
if (ide - 1) % 18 >= 9, calculateExpression1 is used to generate the simualtedExpressionMatrix
if ide < 55, processDataForOptimizing is used to generate the error expression
if 54 < ide < 109, processDataForOptimizingUpdate1 is used to generate the error expression
if ide > 108, processDataForOptimizingBL is used to generate the error expression
all other variations in ide are used to import different .CSV files
number of nodes set in SetMapper can be changed, more nodes will run faster until a max number of nodes required is reach
number of nodes given to SetMapper should equal number of nodes * number of cpus per node which are set in mysticSimulated.sbatch
increasing number of bins given to LatticeSolver will increase accuracy but will also increase time
running a network of 278 genes, 34 TFs, and 44 samples with 4 bins took about 14 hours

mysticSimulatedMain.py -- defines functions transpose, sign, createParams, createParamsUpdate1, createParamsBL, createBounds, calculateExpression, calculatedExpression1, calculatedExpressionBL, and class process with functions processDataForOptimizing, processDataForOptimizingUpdate1, and processDataForOptimizingBL
tranpose takes a list of lists of numbers as input and outputs a tranposed list of lists of numbers
sign takes a list of lists of numbers as inputs and converts all negative numbers to -1, positive numbers to 1, and zeros to 0
createParams, createParamsUpdate1, and createParamsBL take signedBinaryCS and binaryTFA as input and output of list of strings that represent the parameters that are optimized
createBounds takes params (a list of strings) as input and outputs two lists, one for the lower bound and one for the upper bounds
calculateExpression and calculateExpression1 take simulatedControlStrength, simulatedActivity, simulatedBaselines, and simulatedScaling as inputs and output a simulatedExpression list of lists of numbers
calculatedExpressionBL takes simulatedControlStrength and simulatedActivity as input and outputs a simulatedExpression list of lists of numbers
process takes signedBinaryCS, binaryTFA, simulatedExpressionMatrix, and params as input
processDataForOptimizing, processDataForOptimizingUpdate1, and processDataForOptimizingBL all take a list of numbers of length equal to the number of params as input and generate the error between the simulatedExpression matrix and the expressionMatrix calculated by using the input list as parameter values

mysticSimulated.sbatch -- runs mysticSimulated.py, to execute on cluster "sbatch mysticSimulated.sbatch"
set array = to the range of ides that mysticSimulated.py should be run on
value assigned to nodes * value assigned to cpus per node should equal the number of nodes set in SetMapper in mysticSimulated.py

fixParamRulesSimulated.py -- takes in paramRulesSimulated and converts it into a list of rules that Mathematica can use
saves the fixed list of param rules to a file named fixParamRulesSimulated in the folder paramRulesMystic

fixParamRulesSimulated.sbatch run fixParamRulesSimulated.py, to execute on cluster "sbatch fixParamRulesSimulated.sbatch"
