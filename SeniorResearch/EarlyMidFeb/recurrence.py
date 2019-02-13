###############################################################################
# Author : Ian Davidson                                                       #
# Class: CS496 Senior Capstone Project                    	              #
# Date: Feb 10, 2019		  			                      #
# About: This program is intended to run in Python 2.7.15 with librarys numpy,#
#	   collections, and sys.                                              #
# Check notes.txt and ../README.txt	                                      #
###############################################################################


# import sys
# import numpy as np
from collections import *


# dfa encoding from recurrenceRef.pdf hardcoded for simplicity's sake
Transitions = [[2,1],[3,2],[4,2],[1,4],[3,2]]
Start = [1,0,0,0,0]
FINAL = [0,0,0,0,1] #refered to as basecase
inputLength = 100000 #used for later when an iterative solution to this problem Implemented

###############################################################################
# This function recursively executes a recurrence equation                    #
# we reference the global var "Transitions" to generate a recurrence formula  #
# for any state with a given current input left, denoted by "n"               #
# Also note this algorithm is 2^n time complexity,                            #
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
#  alpha = the 2d array used to store values of n * |Q|                       #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
#                                                                             #
###############################################################################
def executeRecurrenceIterative(alpha, n, sigmaSize):
#    print("n: ", n, " sigmaSize", sigmaSize)
    for i in range(1, n):
        for j in range(len(Transitions)):

            alpha[i][j] = float(alpha[i-1][Transitions[j][0]] + alpha[i-1][Transitions[j][1]])/2
            # print("alpha[",i,"][",j,"]: ",  alpha[i][j])
            #print("i: ", i, " j: ", j)
        # ELSE: dont do anything
    return float(float(alpha[n-1][Transitions[0][0]] + alpha[n-1][Transitions[0][1]])/2)






###############################################################################
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
    n = int(input("please enter an integer input size for hardcoded DFA: "))
    alpha = [None] * n
    for i in range(n):
        alpha[i] = [0] * len(Transitions)

    alpha[0] = FINAL
    #print(alpha[0])
    #AcceptingConvergenceProbability = executeRecurrenceRecursive(0, n)
    AcceptingConvergenceProbability = executeRecurrenceIterative(alpha, n, sigmaSize)
    print("Accepting Probability: ", AcceptingConvergenceProbability, " with n=", n)
    #for k in range(len(Transitions)):
        #print("sequence of congergence on state: ", k)
        #for l in range(len(alpha)):
            #print("for input size: ", l, " probability: ", alpha[l][k])
    return
main()
