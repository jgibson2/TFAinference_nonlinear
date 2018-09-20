"""takes a list of lists of numbers as input"""


def transpose(matrix):
    try:
        a = matrix[0][0]
    except:
        return "Error: Incorrect dimensions, need a list of lists of numbers"
    try:
        a = matrix[0][0][0]
        return "Error: Incorrect dimensions, need a list of lists of numbers"
    except:
        pass
    newMatrix = []
    for i in range(len(matrix[0])):
        row = []
        for j in range(len(matrix)):
            row.append(matrix[j][i])
        newMatrix.append(row)
    return (newMatrix)


"""takes a list of lists of numbers as input"""


def sign(matrix):
    try:
        a = matrix[0][0]
    except:
        return "Error: Incorrect dimensions, need a list of lists of numbers"
    try:
        a = matrix[0][0][0]
        return "Error: Incorrect dimensions, need a list of lists of numbers"
    except:
        pass
    newMatrix = []
    for i in range(len(matrix)):
        row = []
        for j in range(len(matrix[0])):
            if matrix[i][j] == 0:
                row.append(0)
            elif matrix[i][j] > 0:
                row.append(1)
            else:
                row.append(-1)
        newMatrix.append(row)
    return newMatrix


"""takes in controlStrength and activity matrices and baseline and scaling vectors and returns an expression matrix using the original non-linear model"""


def calculateExpression(controlStrength, activity, baselines, scaling):
    expressionMatrix = []
    numGenes = len(controlStrength)
    numTFs = len(controlStrength[0])
    numSamples = len(activity[0])
    for i in range(numSamples):
        sample = []
        for j in range(numGenes):
            expressionGeneSample = 0
            """baseline and scaling values are specific to each gene, tfa values are specific to each TF-sample pair, and cs values are specific to each gene-TF pair
            All variable that are being learned on are positive but calculateExpression takes in asigned controlStrength matrix in order to differentiate activators from repressors
            For activators (cs > 0), gene expression = baseline + scaling * tfa / (tfa + cs)
            For repressors (cs < 0), gene expression = baseline + scaling / (tfa + cs)"""
            for k in range(numTFs):
                if controlStrength[j][k] > 0:
                    expressionGeneSample += activity[k][i] / (activity[k][i] + controlStrength[j][k])
                if controlStrength[j][k] < 0:
                    expressionGeneSample += 1 / (activity[k][i] - controlStrength[j][k])
            sample.append(scaling[j] * expressionGeneSample + baselines[j])
        expressionMatrix.append(sample)
    return (transpose(expressionMatrix))


"""takes in controlStrength and activity matrices and baseline and scaling vectors and returns an expression matrix using the updated non-linear model"""


def calculatedExpression1(controlStrength, activity, baselines, scaling):
    expressionMatrix = []
    baseline = []
    numGenes = len(controlStrength)
    numTFs = len(controlStrength[0])
    numSamples = len(activity[0])
    for i in range(numSamples):
        samples = []
        for j in range(numGenes):
            numRepressors = 0
            expressionGeneSample = 0
            """baseline and scaling values are specific to each gene, tfa values are specific to each TF-sample pair, and cs values are specific to each gene-TF pair
            All variables that are being learned on are positive but calculateExpression1 takes in a signed controlStrength matrix in order to differentiate activators from repressors
            For activators (cs > 0), gene expression = scaling * tfa / (tfa + cs)
            For repressors (cs < 0), gene expression = scaling / (tfa +cs)
            If agene does not have an repressors, baseline is added to gene expression
            if a gene has any represssors, no baseline is added"""
            for k in range(numTFs):
                if controlStrength[j][k] < 0:
                    numRepressors += 1
                if controlStrength[j][k] > 0:
                    expressionGeneSample += activity[k][i] / (activity[k][i] + controlStrength[j][k])
                if controlStrength[j][k] < 0:
                    expressionGeneSample += 1 / (activity[k][i] - controlStrength[j][k])
            if numRepressors == 0:
                baseline.append(baselines[j])
            else:
                baseline.append(0)
            samples.append(scaling[j] * expressionGeneSample + baseline[j])
        expressionMatrix.append(samples)
    return (transpose(expressionMatrix))


"""takes in controlStrength and activity matrices and returns a simulated expression matrix using the bi-linear model"""


def calculateExpressionBL(controlStrength, activity):
    expressionMatrix = []
    numGenes = len(controlStrength)
    numTFs = len(controlStrength[0])
    numSamples = len(activity[0])
    for i in range(numSamples):
        sample = []
        for j in range(numGenes):
            expressionGeneSample = 0
            """tfa values are specific to each TF-sample pair and cs values are specific to each gene-TF pair
            All variables that are being learned are positive but calculatedExpressionBL takes in asigned controlStrength matrix in order to differentiate activators from repressors
            For activators (cs > 0), gene expression = cs * tfa
            For repressors (cs < 0), gene expression = -cs * tfa
            Both the controlStrength and activity matrices are given one extra row/column in the TF dimension, so that the extra tfa * cs for each gene-sample pair acts as a baseline"""
            for k in range(numTFs):
                expressionGeneSample += activity[k][i] * controlStrength[j][k]
            sample.append(expressionGeneSample)
        expressionMatrix.append(sample)
    return (transpose(expressionMatrix))


"""returns the list of parameters that are used by the original non-linear model
Need to make it a list of strings because python cannot have a list of symbolic variables"""


def createParams(signedBinaryCS, binaryTFA):
    params = []
    numGenes = len(signedBinaryCS)
    numTFs = len(signedBinaryCS[0])
    numSamples = len(binaryTFA[0])
    for i in range(numSamples):
        for j in range(numGenes):
            if i == 0:
                params.append("s[" + str(j) + "]")
                params.append("b[" + str(j) + "]")
            for k in range(numTFs):
                if j == 0:
                    if binaryTFA[k][i] != 0:
                        params.append("tfa[" + str(k) + "," + str(i) + "]")
                if i == 0:
                    if signedBinaryCS[j][k] != 0:
                        params.append("cs[" + str(j) + "," + str(k) + "]")
    return params


"""returns the list of parameters that are used by the updated non-linear model
Need to make it a list of strings because python cannot have a list of symbolic variables"""


def createParamsUpdate1(signedBinaryCS, binaryTFA):
    params = []
    numGenes = len(signedBinaryCS)
    numTFs = len(signedBinaryCS[0])
    numSamples = len(binaryTFA[0])
    for i in range(numSamples):
        for j in range(numGenes):
            numRepressors = 0
            if i == 1:
                params.append("s[" + str(j) + "]")
            for k in range(numTFs):
                if signedBinaryCS[j][k] < 0:
                    numRepressors += 1
                if j == 0:
                    if binaryTFA[k][i] != 0:
                        params.append("tfa[" + str(k) + "," + str(i) + "]")
                if i == 0:
                    if signedBinaryCS[j][k] != 0:
                        params.append("cs[" + str(j) + "," + str(k) + "]")
            if i == 0:
                if numRepressors == 0:
                    params.append("b[" + str(j) + "]")
    return params


"""returns the list of parameters that are used by the bi-linear model
Need to make it a list of strings because python cannot have a list of symbolic variables"""


def createParamsBL(signedBinaryCS, binaryTFA):
    params = []
    numGenes = len(signedBinaryCS)
    numTFs = len(signedBinaryCS[0])
    numSamples = len(binaryTFA[0])
    for i in range(numSamples):
        for j in range(numGenes):
            for k in range(numTFs):
                if i == 0:
                    params.append("cs[" + str(j) + "," + str(k) + "]")
                if j == 0:
                    params.append("tfa[" + str(k) + "," + str(i) + "]")
    return params


"""creates a list of bounds of length equal to the number of parameters
None makes it so that there is no bound"""


def createBounds(params):
    lb = []
    ub = []
    for param in params:
        ub.append(None)
        if param[0] == "b":
            lb.append(None)
        else:
            lb.append(0.001)
    return (lb, ub)


class process():
    def __init__(self, signedBinaryCS, binaryTFA, geneExpressionMatrix, params):
        self.signedBinaryCS = signedBinaryCS
        self.binaryTFA = binaryTFA
        self.geneExpressionMatrix = geneExpressionMatrix
        self.params = params

    def processDataForOptimizing(self, x):
        errorExpression = 0
        s = {}
        b = {}
        tfa = {}
        cs = {}
        """uses the list of parameters to create dictionaries with each scaling, baseline, tfa, and cs parameters
        Need to do this because python cannot convert the list of parameters from a list of string to a list of variables
        Functions for Mystic to optimize can only have one input, a list of variables, so each dictionary value is assigned to a corresponding Mystic input variable"""
        params = self.params
        signedBinaryCS = self.signedBinaryCS
        binaryTFA = self.binaryTFA
        geneExpressionMatrix = self.geneExpressionMatrix
        for j in range(len(params)):
            i = params[j]
            begin1 = i.index("[") + 1
            end2 = i.index("]")
            if i[0] == "t":
                end1 = i.index(",")
                begin2 = end1 + 1
                tfa[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
            elif i[0] == "c":
                end1 = i.index(",")
                begin2 = end1 + 1
                cs[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
            elif i[0] == "s":
                s[int(i[begin1:end2])] = x[j]
            else:
                b[int(i[begin1:end2])] = x[j]
        numGenes = len(signedBinaryCS)
        numTFs = len(signedBinaryCS[0])
        numSamples = len(binaryTFA[0])
        for i in range(numSamples):
            for j in range(numGenes):
                expressionGeneSample = 0
                """creates an error expression for the optimizer to minimize
                Error expression for a single gene-sample pair equals the difference between the actual expression value and the expression equation of variables all squared
                The total error expression that the optimizer uses is the sum of the error expression over all genes and samples"""
                for k in range(numTFs):
                    if binaryTFA[k][i] > 0 and signedBinaryCS[j][k] > 0:
                        expressionGeneSample += tfa[k, i] / (tfa[k, i] + cs[j, k])
                    elif binaryTFA[k][i] > 0 and signedBinaryCS[j][k] < 0:
                        expressionGeneSample += 1 / (tfa[k, i] + cs[j, k])
                    if binaryTFA[k][i] == 0 and signedBinaryCS[j][k] < 0:
                        expressionGeneSample += 1 / cs[j, k]
                expressionGeneSample *= s[j]
                expressionGeneSample += b[j]
                errorExpression += (geneExpressionMatrix[j][i] - expressionGeneSample) ** 2
        return (errorExpression)

    """creates the error expression with the updated non-linear model"""

    def processDataForOptimizingUpdate1(self, x):
        params = self.params
        signedBinaryCS = self.signedBinaryCS
        binaryTFA = self.binaryTFA
        geneExpressionMatrix = self.geneExpressionMatrix
        errorExpression = 0
        s = {}
        b = {}
        tfa = {}
        cs = {}
        """uses the list of parameters to create dictionaries with each scaling, baseline, tfa, and cs parameter
        Need to do this because python cannot convert the list of parameters from a list of string to a list of variables
        Functions for Mystic to optimize can only have on input, a list of variables, so each dictionary value is assigned to a corresponding Mystic input variable"""
        for j in range(len(params)):
            i = params[j]
            begin1 = i.index("[") + 1
            end2 = i.index("]")
            if i[0] == "t":
                end1 = i.index(",")
                begin2 = end1 + 1
                tfa[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
            elif i[0] == "c":
                end1 = i.index(",")
                begin2 = end1 + 1
                cs[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
            elif i[0] == "s":
                s[int(i[begin1:end2])] = x[j]
            else:
                b[int(i[begin1:end2])] = x[j]
        numGenes = len(signedBinaryCS)
        numTFs = len(signedBinaryCS[0])
        numSamples = len(binaryTFA[0])
        for i in range(numSamples):
            for j in range(numGenes):
                numRepressors = 0
                expressionGeneSample = 0
                """creates an error expression for the optimizer to minimize
                Error expression for asingle gene-sample pair equals the difference between the actual expression value and the expression equation of variables all squared
                The total error expression that the optimizer uses is the sum of the error expression over all genes and samples"""
                for k in range(numTFs):
                    if signedBinaryCS[j][k] < 0:
                        numRepressors += 1
                    if binaryTFA[k][i] > 0 and signedBinaryCS[j][k] > 0:
                        expressionGeneSample += tfa[k, i] / (tfa[k, i] + cs[j, k])
                    elif binaryTFA[k][i] > 0 and signedBinaryCS[j][k] < 0:
                        expressionGeneSample += 1 / (tfa[k, i] + cs[j, k])
                    elif binaryTFA[k][i] == 0 and signedBinaryCS[j][k] < 0:
                        expressionGeneSample += 1 / cs[j, k]
                expressionGeneSample *= s[j]
                if numRepressors == 0:
                    expressionGeneSample += b[j]
                errorExpression += (geneExpressionMatrix[j][i] - expressionGeneSample) ** 2
        return errorExpression

    """creates the error expression with the bi-linear model"""

    def processDataForOptimizingBL(self, x):
        errorExpression = 0
        tfa = {}
        cs = {}
        """uses the list of parameters to create dictionaries with each tfa and cs parameter
        Need to do this because python cannot convert the list of parameters from a list of strings to a list of variabels
        Function for Mystic to optimize can only have one input, a list of variables, so each dictionary value is assigned to a corresponding Mystic input variable"""
        for j in range(len(self.params)):
            i = self.params[j]
            begin1 = i.index("[") + 1
            end2 = i.index("]")
            if i[0] == "t":
                end1 = i.index(",")
                begin2 = end1 + 1
                tfa[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
            elif i[0] == "c":
                end1 = i.index(",")
                begin2 = end1 + 1
                cs[int(i[begin1:end1]), int(i[begin2:end2])] = x[j]
        numGenes = len(self.signedBinaryCS)
        numTFs = len(self.signedBinaryCS[0])
        numSamples = len(self.binaryTFA[0])
        for i in range(numSamples):
            for j in range(numGenes):
                expressionGeneSample = 0
                """creates an error expression for the optimizer to minimize
                Error expression for a single gene-sample pair equals the difference between the actual expression value and the expression equation of variables all squared
                The total error expression that the optimizer uses is the sum of the error expression over all genes and samples"""
                for k in range(numTFs):
                    if self.signedBinaryCS[j][k] != 0:
                        if self.signedBinaryCS[j][k] > 0:
                            expressionGeneSample += cs[j, k] * tfa[k, i]
                        else:
                            expressionGeneSample -= cs[j, k] * tfa[k, i]
                errorExpression += (self.geneExpressionMatrix[j][i] - expressionGeneSample) ** 2
        return (errorExpression)
