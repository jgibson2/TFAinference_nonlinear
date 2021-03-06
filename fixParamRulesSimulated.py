"""ide is only used for accessing the paramRules of different simulated data"""
import sys

ide = int(sys.argv[1])

"""because python does not have a rules built in, the output of the optimizer has to be as list of strings
Also, python lists are denoted with "[" and "]", but Mathematica lists are denotes with "{" and "}", so this converts the list of strings into a single string that can be imported in Mathematica and casted to a list of rules by replacing the opening and closing "[" and "]" of the list with opening and closing "{" and "}" and removing linebreaks
The string is then save to a file that can be imported into Mathematica"""
try:
    with open("paramRulesMystic//paramRulesSimulated" + str(ide)) as f:
        lst = f.read()
    newLst = "{"
    for i in range(len(lst)):
        j = lst[i]
        if j != "\'" and j != "[" and j != "]":
            newLst += str(j)
        if (j == "[" or j == "]") and i != 0 and i != len(lst) - 1:
            newLst += str(j)
    newLst += "}"
except:
    newLst = "0"
with open("paramRulesMystic//fixParamRulesSimulated" + str(ide), "w") as file:
    file.write(newLst)
