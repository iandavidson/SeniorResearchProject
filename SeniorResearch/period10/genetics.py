

import numpy as np
import sys
from collections import *
import copy as copy
import math

# StateTitles = ["AA", "Aa", "aa"]
# Transitions = [[1,0,0],[0.25, .5, .25],[0,0,1]]
# stateSize = 3


StateTitles = ["AA, AA","aa, aa","AA, aA","aa, aA","aA, aA","aa, AA"]
Transitions = [
            [1,0,0,0,0,0], #AA, AA
            [0,1,0,0,0,0], #aa, aa
            [0.25,0,0.5,0,0,0.25], #AA, aA
            [0,0.25,0,0.5,0,0.25], #aa, aA
            [0,0,0,0,0,1], #aA, aA
            [0.0625,0.0625,0.25,0.25,0.125,0.25]] #aa, AA
            #these transitions represent the sibmating markov chain,
            #basically we are saying we make a punitt square with the state labels, then we randomly choose two results with repreition
stateSize = 6
wordLen = 50 # in turn this means 50 generations


def ComputeDensityMatrix(Final, start):
    q = len(Transitions)
    #create square matrix of transition table
    delta = copy.deepcopy(Transitions)
    # print delta
    #use the A^2 trick. use the odd and even cases with recursive function.


    newDelta = divideByConstantHelper(delta, wordLen)
    # print "new delta:"
    # print newDelta[start]

    #do: dotProduct(newDelta.FirstRow, FINAL) <- assuming the first state is starting
    acceptedSum = np.dot(newDelta[start], Final)
    return acceptedSum

#  does log(n) arithmatic steps in computing (delta)^(n)
def divideByConstantHelper(delta, n):
    if n == 0:
        return np.identity(len(delta))#identity matrix

    elif n % 2 == 0: #n is even
        newDelta = copy.deepcopy(divideByConstantHelper(delta, n/2))
        return np.matmul(newDelta, newDelta)
    else: # n % 2 == 1
        newDelta = copy.deepcopy(divideByConstantHelper(delta, (n-1)/2))
        return np.matmul(np.matmul(newDelta, newDelta), delta)


def ComputeDensityRecurrence(Final, Start):
#assuming there is only going to be one accepting state.
    old = copy.deepcopy(Final)
    current = [None] * stateSize
    for i in range(wordLen + 1):
        for j in range(stateSize):

            Sum = 0.0

            for k in range(stateSize):
                #iterate through all of the transitions of

                Sum += Transitions[j][k] * old[k]

                # Sum += Transitions[j][k] * old[j]
            #assign current[j] from sum

            current[j] = Sum
        old = copy.deepcopy(current)

    # return statement looks like this because start state is 0
    # return float(old[Transitions[Start][0]] + old[Transitions[Start][1]])/2.0

    return current[Start]



def main():
    Final = [0] * stateSize
    # Start = [0] * 6
    for i in range(stateSize):
        print "Computing convergence over " + str(wordLen) + " generations for 2 children with the alleles: " + str(StateTitles[i])
        Final[i] = 1
        for j in range(stateSize):
            # print Final
            print "\t Assuming we start with parent traits : " + str(StateTitles[j])
            result = ComputeDensityMatrix(Final, j)
            # ConvergencePercent = ComputeDensityRecurrence(Final, j)
            print "\t Convergence probability is: " + str(result)
            print "\n"
        Final[i] = 0
        print "\n"
    # Final = [0,1,0]
    # Start = 1

    # print ComputeDensityRecurrence(Final, Start)

    # Final = [0,1,0]


main()
