import numpy as np
import sys
from collections import defaultdict
import copy as copy
import matplotlib.pyplot as plt
import time as time


###############################################################################
# MarkovChain class: to represent a regular language for usage for the data model : markovchain
# members: 
#     States -> number of states in instance                            
#     Transitions -> uses a 2d list of States by States dimentions. #
#
# Operations of class:
#   Operations to compute density for single value:
#     ComputeDensityMatrixSingleStart // this function requires a singlar start state (100% starting in a single state)
#     ComputeDensityRecurrenceSingleStart // this function requires a singlar start state (100% starting in a single state)
#     ComputeDensityMatrixVectorStart //  this function requires a vector-start array (distributed starting state)
#     ComputeDensityRecurrenceVectorStart // this function requires a vector-start array (distributed starting state)
#
#   Operations to compute density for a range of values:
#     ComputeDensityRangeRecurrence // this function requires a vector-start array (distributed starting state)
# 
# 
###############################################################################
class MarkovChain:
    def __init__(self, states, transitions):
        self.States = states
        self.Transitions = transitions


    def updateInstance(self, states,transitions):
        self.States  = states
        self.Transitions = transitions


    def ComputeDensityMatrixVectorStart(self, Final, Start, wordLen):
#     Start -> list with States elements, where the sum => 1, provides initial conditions.
#     Final -> list with States elements, where values are 0 || 1 -> for states with indexs we are interested in convergence on. #
#               (same as final vector for DFA/NFA)          
#     wordLen -> length of word we are computing density for.
#     Return: the convergence/density at wordLen input size
        q = len(self.Transitions)
        #create square matrix of transition table
        delta = copy.deepcopy(self.Transitions)
        # print delta
        #use the A^2 trick. use the odd and even cases with recursive function.


        newDelta = self.divideByConstantHelper(delta, wordLen)
        # print "new delta:"
        # print newDelta[start]

        #newDelta = Start * newDelta  -> here we do matrix mult for [1 x n]*[n x n] => [1 x n]
        newDelta = np.matmul(Start, newDelta)  

        #acceptedSum = dot product of ([1 x n], [n x 1]) -> result is a number between 0 and 1.
        acceptedSum = np.dot(newDelta, Final)

        return acceptedSum




    def ComputeDensityMatrixSingleStart(self, Final, start, wordLen):
#     Start -> is an int of the state index for 100% starting state.
#     Final -> list with States elements, where values are 0 || 1 -> for states with indexs we are interested in convergence on. #
#               (same as final vector for DFA/NFA)          
#     wordLen -> length of word we are computing density for.
#     Return: the convergence/density at wordLen input size
        q = len(self.Transitions)
        #create square matrix of transition table
        delta = copy.deepcopy(self.Transitions)
        # print delta
        #use the A^2 trick. use the odd and even cases with recursive function.


        newDelta = self.divideByConstantHelper(delta, wordLen)



        #do: dotProduct(newDelta.FirstRow, FINAL) <- assuming <start>  is starting index
        acceptedSum = np.dot(newDelta[start], Final)

        return acceptedSum

    #  does log(n) arithmatic steps in computing (delta)^(n)
    def divideByConstantHelper(self, delta, n):
        if n == 0:
            return np.identity(len(delta))#identity matrix

        elif n % 2 == 0: #n is even
            newDelta = copy.deepcopy(self.divideByConstantHelper(delta, n/2))
            return np.matmul(newDelta, newDelta)
        else: # n % 2 == 1
            newDelta = copy.deepcopy(self.divideByConstantHelper(delta, (n-1)/2))
            return np.matmul(np.matmul(newDelta, newDelta), delta)


    def ComputeDensityRecurrenceSingleStart(self, Final, Start, wordLen):
#     Start -> Start -> is an int of the state index for 100% starting state.
#     Final -> list with States elements, where values are 0 || 1 -> for states with indexs we are interested in convergence on. #
#               (same as final vector for DFA/NFA)          
#     wordLen -> length of word we are computing density for.
#     Return: the convergence/density at wordLen input size
        old = copy.deepcopy(Final)
        current = [None] * self.States
        for i in range(wordLen + 1):
            for j in range(self.States):

                Sum = 0.0

                for k in range(self.States):
                    #iterate through all of the transitions of

                    Sum += self.Transitions[j][k] * old[k]

                    # Sum += Transitions[j][k] * old[j]
                #assign current[j] from sum

                current[j] = Sum
            old = copy.deepcopy(current)

        return current[Start]


    def ComputeDensityRecurrenceVectorStart(self, Final, Start, wordLen):
#     Start -> list with States elements, where the sum => 1, provides initial conditions.
#     Final -> list with States elements, where values are 0 || 1 -> for states with indexs we are interested in convergence on. #
#               (same as final vector for DFA/NFA)          
#     wordLen -> length of word we are computing density for.
#     Return: the convergence/density at wordLen input size
        old = copy.deepcopy(Final)
        current = [None] * self.States
        for i in range(wordLen):
            for j in range(self.States):

                Sum = 0.0

                for k in range(self.States):
                    #iterate through all of the transitions of each state

                    Sum += self.Transitions[j][k] * old[k]

                #assign current[j] from sum

                current[j] = Sum
            old = copy.deepcopy(current)

        convergence = 0.0
        for i in range(self.States):
            convergence += (current[i] * Start[i])
        return convergence


    def ComputeDensityRangeRecurrence(self, Final, Start, idxStart, idxFinish):
#     Start -> list with States elements, where the sum => 1, provides initial conditions.
#     Final -> list with States elements, where values are 0 || 1 -> for states with indexs we are interested in convergence on. #
#               (same as final vector for DFA/NFA)          
#     idxStart -> the input length at which to start logging the convergence Based off Start and Finish vectors
#     idxFinish -> the input length at which we stop logging the convergence 
#     Return: an array of convergence/density from input size idxStart to idxFinish
        convergence = 0.0
        log = []
        old = copy.deepcopy(Final)
        current = [None] * self.States
        for i in range(idxFinish):
            for j in range(self.States):

                Sum = 0.0

                for k in range(self.States):
                    #iterate through all of the transitions of each state

                    Sum += self.Transitions[j][k] * old[k]

                    # Sum += Transitions[j][k] * old[j]
                #assign current[j] from sum

                current[j] = Sum
            old = copy.deepcopy(current)
            
            if i >= idxStart:
                convergence = 0.0
                for i in range(self.States):
                    convergence += (current[i] * Start[i])
                log.append( convergence)


        convergence = 0.0
        for i in range(self.States):
            convergence += (current[i] * Start[i])
        log.append(convergence)
        return log

def main():


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
    markovChain = MarkovChain(stateSize, Transitions)

    final = [1,1,0,0,0,0]
    start = 5


    print "ComputeDensityMatrixSingleStart:"
    ComputeDensityMatrixSingleStartResult = markovChain.ComputeDensityMatrixSingleStart(final, start, 6)
    print str(ComputeDensityMatrixSingleStartResult) + "\n"

    
    print "ComputeDensityRecurrenceSingleStart:"
    ComputeDensityRecurrenceSingleStartResult = markovChain.ComputeDensityRecurrenceSingleStart(final, start, 5)
    print str(ComputeDensityRecurrenceSingleStartResult) + "\n"


    startVector = [0.25, 0.25, 0, 0, .25, .25]


    print "ComputeDensityMatrixVectorStart:"
    ComputeDensityMatrixVectorStartResult = markovChain.ComputeDensityMatrixVectorStart(final, startVector, 5)
    print str(ComputeDensityMatrixVectorStartResult) + "\n"


    print "ComputeDensityRecurrenceVectorStart:"
    ComputeDensityRecurrenceVectorStartResult = markovChain.ComputeDensityRecurrenceVectorStart(final, startVector, 5)
    print str(ComputeDensityRecurrenceVectorStartResult) + "\n"

    print "ComputeDensityRangeRecurrence: with start: 1, finish: 5" 
    ComputeDensityRangeRecurrenceResult = markovChain.ComputeDensityRangeRecurrence(final, startVector, 1, 5)
    print ComputeDensityRangeRecurrenceResult

main()  

# ComputeDensityRangeRecurrence // this function requires a vector-start array (distributed starting state)




# StateTitles = ["AA, AA","aa, aa","AA, aA","aa, aA","aA, aA","aa, AA"]
# Transitions = [
#             [1,0,0,0,0,0], #AA, AA
#             [0,1,0,0,0,0], #aa, aa
#             [0.25,0,0.5,0,0,0.25], #AA, aA
#             [0,0.25,0,0.5,0,0.25], #aa, aA
#             [0,0,0,0,0,1], #aA, aA
#             [0.0625,0.0625,0.25,0.25,0.125,0.25]] #aa, AA
#             #these transitions represent the sibmating markov chain,
#             #basically we are saying we make a punitt square with the state labels, then we randomly choose two results with repreition
# stateSize = 6
# wordLen = 50 