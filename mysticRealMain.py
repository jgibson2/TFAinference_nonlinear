"""takes a list of lists of numbers as input"""


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
