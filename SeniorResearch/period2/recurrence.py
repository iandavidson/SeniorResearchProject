###############################################################################
# Author : Ian Davidson                                                       #
# Class: CS496 Senior Capstone Project                    	                  #
# Date: Feb 16, 2019		  			                                      #
# About: This program is intended to run in Python 2.7.15 with librarys numpy,#
#	   collections, and sys.                                                  #
# Check notes.txt and ../README.txt	                                          #
###############################################################################

import numpy as np
import sys
from collections import *
import copy

# dfa encoding from recurrenceRef.pdf hardcoded for simplicity's sake
Transitions = [[1,5],[3,5],[0,5],[4,5],[4,3],[6,6],[5,6]]
Start = [1,0,0,0,0,0,0]
FINAL = [0,0,0,0,0,0,1]
inputLength = 100050 #used for later when an iterative solution to this problem Implemented




def divideByConstantFactor(n):
    q = len(Transitions)
    #create square matrix of transition table
    delta = np.zeros((q,q))

    #make square matrix from transitions
    for j in range(q):
        for k in range(2):
            delta[j][Transitions[j][k]] += 1

    #use the A^2 trick. use the odd and even cases with recursive function.
    newDelta = divideByConstantHelper(delta, n)


    #do: start * newDelta * FINAL
    acceptedSum = np.matmul(np.matmul(Start, newDelta), FINAL)

    #then we need to compute 2^n
    denominator = 2**n

    #then do acceptedSum/2^n then return product
    # print "accepted # of strings: " + str(acceptedSum)
    # print "total possible string: " + "2^" + str(n)

    return float(acceptedSum/denominator)



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







###############################################################################
# This function recursively executes a recurrence equation                    #
# we reference the global var "Transitions" to generate a recurrence formula  #
# for any state with a given current input left, denoted by "n"               #
# Also note this algorithm is 2^n time complexity,                            #
# This is known as Version1 for computing recurrence based off transitions.   #
###############################################################################
def executeRecurrenceRecursive(i, n):
    if(n == 0): #basecase
        return FINAL[i]

    temp = 0.0
    for k in range(len(Transitions[i])): #tells us how many transitions are there out of state i to get its respective relation
        probability = float(executeRecurrenceRecursive(Transitions[i][k], n-1))
        #print("prob for i= ", i, " n = ", n, " : ", probability)
        temp += probability
    return float((temp)/len(Transitions[i]))



###############################################################################
#  alpha = the 2d array used to store values of n * |Q|   (very bad space complexity-wise)#
#  |Q| => len(Transitions[anyvalue: 0->n])                                    #
#  n = the value of F[0](n) is what we are computing                          #
#  sigmaSize => len(Transitions[0]) - set size of input alphabet              #
#  This is version2                                                           #
#                                                                             #
###############################################################################
def executeRecurrenceIterative(alpha, n, sigmaSize):
#    print("n: ", n, " sigmaSize", sigmaSize)
    for i in range(1, n):
        for j in range(len(Transitions)):

            alpha[i][j] = float(alpha[i-1][Transitions[j][0]] + alpha[i-1][Transitions[j][1]])/2
            print("alpha[",i,"][",j,"]: ",  alpha[i][j])
            #print("i: ", i, " j: ", j)
        # ELSE: dont do anything
    return float(float(alpha[n-1][Transitions[0][0]] + alpha[n-1][Transitions[0][1]])/2)


###############################################################################
#  |Q| => len(Transitions[anyvalue: 0->|Q|])                                    #
#  n = the value of F[0](n) is what we are computing                          #
#  sigmaSize => len(Transitions[0]) - set size of input alphabet              #
#  Version3 Key difference from V2 are, better space complexity -> O(sigmaSize * |States|)                                                                           #
#                                                                             #
###############################################################################
def executeRecurrenceIterativeVersion2(n, sigmaSize):

    old = copy.deepcopy(FINAL)
    current = [None] * len(Transitions)

    for i in range(n):
        for j in range(len(Transitions)):
            # print "old[Transitions[" +str(j) + "][0]]: " + str(old[Transitions[j][0]]) + ";  old[Transitions[" +str(j) + "][1]]:" + str(old[Transitions[j][1]])
            current[j] = float(float(old[Transitions[j][0]] + old[Transitions[j][1]]))/2

            if FINAL[j] == 1 and i > 100000:
                print "state: "+ str(j) + "; n: " + str(i) + "; " + str(current[j])
            # print "old[Transitions[" + str(j) + "][0]]: ", old[Transitions[j][0]]
            # print "old[Transitions[" + str(j) + "][1]]: ", old[Transitions[j][1]]
            #print("i: ", i, " j: ", j)
        # ELSE: dont do anything
        # print "iteration: ", i
        # print "old: ", old
        # print "current: ", current
        old = copy.deepcopy(current)

    # return statement looks like this because start state is 0
    return float(float(old[Transitions[0][0]] + old[Transitions[0][1]])/2)




######### ######################################################################
# This program will do 2 steps with a hardcoded encoding of a DFA:            #
#	  Step1:   replace all the transitions by transitions of a weight = 1/2   #
#	          This represents the dfa as a Partially Observable Markov Chain  #
#     Step2:   write equations of the form; each state gets an equation       #
#	          P[i](n) = SUM(P[j](n-1)) | j = delta(i,a) for some a;           #	                                                                          #
#	 Currently there my function executeRecurrenceRecursive() which uses      #
#   DFA encoding stored, to generate and execute respective recurrence formula#
###############################################################################
def main():
    sigmaSize = len(Transitions[0]) # is size of sigma

    # len(Transitions) is the number of states
    n =  int(input("please enter an integer input size for hardcoded DFA: "))
    if(type(n) is not int):
        exit(1)


    print "Matrix method density at    n: " + str(divideByConstantFactor(n))



    AcceptingConvergenceProbability = executeRecurrenceIterativeVersion2(n, sigmaSize)
    print "Iterative method density at n: " + str(AcceptingConvergenceProbability)



    # AcceptingConvergenceProbability = executeRecurrenceIterative(alpha, n, sigmaSize)
    # print("Accepting Probability V1: ", AcceptingConvergenceProbability, " with n=", n)

    # for k in range(len(Transitions)):
    #     if(FINAL[k] == 1):
    #         print "sequence of congergence on state: ", k
    #         for l in range(len(alpha)):
    #             print "for input size: ", l, " probability: ", alpha[l][k]


    return



main()
