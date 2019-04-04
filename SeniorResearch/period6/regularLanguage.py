###############################################################################
# Author : Ian Davidson                                                       #
# Class: CS496 Senior Capstone Project                    	                  #
# Date: March 5, 2019		  			                                      #
# About: This program is intended to run in Python 2.7.15 with librarys numpy,#
#	   collections, and sys.                                                  #
# Check notes.txt and ../README.txt	                                          #
###############################################################################

#Goal of this file:
    #1 Start with an instance of a dfa.
        #i. Find it's strongly connected components done.
    #2 Create adjacency list for each subchain induced by each SCC/terminal component.
        #i. I was thinking of making a class for each subchain or make subchains members of DFA class.
    #3 Use subchains to run through algorithms Compute-Limit and Compute-Period.


import numpy as np
import sys
from collections import defaultdict
import copy as copy
import matplotlib.pyplot as plt
import time as time

###############################################################################
# DFA class: to represent a regular language for usage in the                 #
# algorithm compute-period                                                    #
# other alogithms implemented for the class are bfs and finding SCCs   #
# PreProcessDFA eliminates useless states. Is explicitly called on an instance#
###############################################################################
class NFA:

    def __init__(self, states, inputSize, final, transitions=None):
        self.States = states
        if transitions == None:
            self.Transitions = defaultdict(list)
            for i in range(states):
                self.Transitions[i] = []
                for j in range(inputSize):
                    self.Transitions[i].append([])
        else:
            self.Transitions = transitions
        self.InputSize = inputSize
        self.Final = final
        self.StronglyCCs = []

    #function used to add an edge directed from stateIdx: <fromState> -> <endState>
    def addEdge(self, fromState, endState, edgeValue):

        self.Transitions[fromState][edgeValue].append(endState)



    # We want to see if there are any unreachable states from the start
    # and remove them, we check at multiple
    def PreProcessDFA(self):

        #call bfs from starting state (always asumming 0 is starting state)
        visited = self.InitialBreadthFirstSearch(0)

        #transform self.Transitions in some way
        correctedTransitions = copy.deepcopy(self.Transitions)
        correctedStates = self.States
        correctedFinal = []
        #remove all the states that are not reachable
        for i in range(self.States):
            if visited[i] == 0:
                #remove state even if it is accepting, dont care because we if unreachable
                correctedTransitions.pop(i)

                correctedStates -= 1
            else:
                correctedFinal.append(self.Final[i])



        if correctedStates == self.States:
            return #no adjustment was required, didn't remove any states

        #assign new adjusted class members back onto self instance

        self.Transitions = correctedTransitions
        self.States = correctedStates
        self.Final = correctedFinal

        return

    def delta(self, state, inputIndex):

        #will return a list!
        return self.Transitions[state][inputIndex]


    def InitialBreadthFirstSearch(self, firstState):

        # Mark all the vertices as not visited
        visited = [0] * (self.States)

        # Create a queue for BFS
        queue = []

        # Mark the source node as
        # visited and enqueue it
        queue.append(firstState)
        visited[firstState] = True

        while queue:

            # Dequeue a vertex from
            # queue and print it
            state = queue.pop(0)

            # Get all adjacent vertices of the
            # dequeued vertex s. If a adjacent
            # has not been visited, then mark it
            # visited and enqueue it
            for i in (self.Transitions[state]):
                for j in range(self.Transitions[state][i]):
                    if visited[self.Transitions[state][i][j]] == 0:
                        queue.append(self.Transitions[state][i][j])
                        visited[self.Transitions[state][i][j]] = 1

        return visited


    def LogDensityRecurrence(self, wordLen, startIndex, endIndex):
        old = copy.deepcopy(self.Final)
        current = [None] * len(self.Transitions)
        log = []
        sum = 0.0
        for i in range(wordLen):
            for j in range(len(self.Transitions)):
                # print "old[Transitions[" +str(j) + "][0]]: " + str(old[Transitions[j][0]]) + ";  old[Transitions[" +str(j) + "][1]]:" + str(old[Transitions[j][1]])

                input0 = 0.0
                delta0 = self.delta(j, 0)
                for k in range(len(delta0)):
                    input0 = float(input0 + old[delta0[k]])
                input0 /= len(delta0)
                input1 = 0.0
                delta1 = self.delta(j, 1)
                for l in range(len(delta1)):
                    input1 = float(input1 + old[delta1[l]])
                input1 /= len(delta1)
                current[j] = float(input0 + input1)/2.0 #force a cast to float type

                # if self.Final[j] == 1 and i > 100000:
                if j == 0 and i <= endIndex and i >= startIndex:
                    #print "state: "+ str(j) + "; n: " + str(i) + "; " + str(current[j])
                    sum += current[j]

            old = copy.deepcopy(current)
            if i <= endIndex and i >= startIndex:
                log.append(sum)
                sum = 0.0

        return log

    def ComputeDensityRecurrence(self, wordLen):
        old = copy.deepcopy(self.Final)
        current = [None] * len(self.Transitions)
        for i in range(wordLen):
            for j in range(len(self.Transitions)):
                # print "old[Transitions[" +str(j) + "][0]]: " + str(old[Transitions[j][0]]) + ";  old[Transitions[" +str(j) + "][1]]:" + str(old[Transitions[j][1]])
                current[j] = float(old[self.Transitions[j][0]] + old[self.Transitions[j][1]])/2.0

            old = copy.deepcopy(current)

        # return statement looks like this because start state is 0
        return float(old[self.Transitions[0][0]] + old[self.Transitions[0][1]])/2.0

    def ComputeDensityMatrix(self, wordLen):
        q = len(self.Transitions)
        #create square matrix of transition table
        delta = np.zeros((q,q))

        #make square matrix from transitions
        for j in range(q):
            for k in range(2):
                for l in range(len(self.Transitions[j][k])): #could even be 0
                    delta[j][self.Transitions[j][k]] += 1

        #use the A^2 trick. use the odd and even cases with recursive function.
        newDelta = self.divideByConstantHelper(delta, wordLen)


        #do: dotProduct(newDelta.FirstRow, FINAL) <- assuming the first state is starting
        acceptedSum = np.dot(newDelta[0], self.Final)

        #then we need to compute 2^n
        denominator = 2**wordLen

        #then do acceptedSum/2^n then return product
        # print "accepted # of strings: " + str(acceptedSum)
        # print "total possible string: " + "2^" + str(n)

        return float(acceptedSum/denominator)



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







###############################################################################
# DFA class: to represent a regular language for usage in the                 #
# algorithm compute-period                                                    #
# other alogithms implemented for the class are bfs and finding SCCs   #
# PreProcessDFA eliminates useless states. Is explicitly called on an instance#
###############################################################################
class DFA:

    def __init__(self, states, inputSize, final, transitions=None):
        self.States = states
        if transitions == None:
            self.Transitions = defaultdict(list)
        else:
            self.Transitions = transitions
        self.InputSize = inputSize
        self.Final = final
        self.StronglyCCs = []

    #function used to add an edge directed from stateIdx: <fromState> -> <endState>
    def addEdge(self, fromState, endState):
        if(len(self.Transitions[fromState]) < self.InputSize):
            self.Transitions[fromState].append(endState)
        else:
            self.Transitions[fromState].pop()
            self.Transitions[fromState].append(endState)
            #either do self or just dont add


    # We want to see if there are any unreachable states from the start
    # and remove them, we check at multiple
    def PreProcessDFA(self):

        #call bfs from starting state (always asumming 0 is starting state)
        visited = self.InitialBreadthFirstSearch(0)

        #transform self.Transitions in some way
        correctedTransitions = copy.deepcopy(self.Transitions)
        correctedStates = self.States
        correctedFinal = []
        #remove all the states that are not reachable
        for i in range(self.States):
            if visited[i] == 0:
                #remove state even if it is accepting, dont care because we if unreachable
                correctedTransitions.pop(i)

                correctedStates -= 1
            else:
                correctedFinal.append(self.Final[i])



        if correctedStates == self.States:
            return #no adjustment was required, didn't remove any states

        #assign new adjusted class members back onto self instance

        self.Transitions = correctedTransitions
        self.States = correctedStates
        self.Final = correctedFinal

        return

    def delta(self, state, inputIndex):

        # if(self.Transitions[state][inputIndex] == null):
            #generate a transition and assign on that ^^
                #generate with: (random-seed * state + inputIndex) % self.states
        return self.Transitions[state][inputIndex]





    def InitialBreadthFirstSearch(self, firstState):

        # Mark all the vertices as not visited
        visited = [0] * (self.States)

        # Create a queue for BFS
        queue = []

        # Mark the source node as
        # visited and enqueue it
        queue.append(firstState)
        visited[firstState] = True

        while queue:

            # Dequeue a vertex from
            # queue and print it
            state = queue.pop(0)

            # Get all adjacent vertices of the
            # dequeued vertex s. If a adjacent
            # has not been visited, then mark it
            # visited and enqueue it
            for i in self.Transitions[state]:
                if visited[i] == 0:
                    queue.append(i)
                    visited[i] = 1

        return visited


    def LogDensityRecurrence(self, wordLen, startIndex, endIndex):
        old = copy.deepcopy(self.Final)
        current = [None] * len(self.Transitions)
        log = []
        sum = 0.0
        for i in range(wordLen):
            for j in range(len(self.Transitions)):
                # print "old[Transitions[" +str(j) + "][0]]: " + str(old[Transitions[j][0]]) + ";  old[Transitions[" +str(j) + "][1]]:" + str(old[Transitions[j][1]])
                current[j] = float(old[self.Transitions[j][0]] + old[self.Transitions[j][1]])/2.0

                # if self.Final[j] == 1 and i > 100000:
                if j == 0 and i <= endIndex and i >= startIndex:
                    #print "state: "+ str(j) + "; n: " + str(i) + "; " + str(current[j])
                    sum += current[j]
                # print "old[Transitions[" + str(j) + "][0]]: ", old[Transitions[j][0]]
                # print "old[Transitions[" + str(j) + "][1]]: ", old[Transitions[j][1]]
                #print("i: ", i, " j: ", j)
            # ELSE: dont do anything
            # print "iteration: ", i
            # print "old: ", old
            # print "current: ", current
            old = copy.deepcopy(current)
            if i <= endIndex and i >= startIndex:
                log.append(sum)
                sum = 0.0
        # return statement looks like this because start state is 0
        # log.append(float(old[self.Transitions[0][0]] + old[self.Transitions[0][1]])/2.0)
        return log

    def ComputeDensityRecurrence(self, wordLen):
        old = copy.deepcopy(self.Final)
        current = [None] * len(self.Transitions)
        for i in range(wordLen):
            for j in range(len(self.Transitions)):
                # print "old[Transitions[" +str(j) + "][0]]: " + str(old[Transitions[j][0]]) + ";  old[Transitions[" +str(j) + "][1]]:" + str(old[Transitions[j][1]])
                current[j] = float(old[self.Transitions[j][0]] + old[self.Transitions[j][1]])/2.0

            old = copy.deepcopy(current)

        # return statement looks like this because start state is 0
        return float(old[self.Transitions[0][0]] + old[self.Transitions[0][1]])/2.0

    def ComputeDensityMatrix(self, wordLen):
        q = len(self.Transitions)
        #create square matrix of transition table
        delta = np.zeros((q,q))

        #make square matrix from transitions
        for j in range(q):
            for k in range(2):
                delta[j][self.Transitions[j][k]] += 1

        #use the A^2 trick. use the odd and even cases with recursive function.
        newDelta = self.divideByConstantHelper(delta, wordLen)


        #do: dotProduct(newDelta.FirstRow, FINAL) <- assuming the first state is starting
        acceptedSum = np.dot(newDelta[0], self.Final)

        #then we need to compute 2^n
        denominator = 2**wordLen

        #then do acceptedSum/2^n then return product
        # print "accepted # of strings: " + str(acceptedSum)
        # print "total possible string: " + "2^" + str(n)

        return float(acceptedSum/denominator)



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






    def findSCCs(self):
        #unimplemented
        return

# 1.  For each vertex u of the graph, mark u as unvisited. Let L be empty.
# 2.  For each vertex u of the graph do Visit(u), where Visit(u) is the recursive subroutine:
#         If u is unvisited then:
#             1.Mark u as visited.
#             2.For each out-neighbour v of u, do Visit(v).
#             3.Prepend u to L.
#         Otherwise do nothing.
# 3.  For each element u of L in order, do Assign(u,u) where Assign(u,root) is the recursive subroutine:
#         If u has not been assigned to a component then:
#             1.Assign u as belonging to the component whose root is root.
#             2.For each in-neighbour v of u, do Assign(v,root).
#         Otherwise do nothing.


def callWithDFAExample0():
    states = 5
    inputSize = 2
    final = [0, 0, 0, 0, 1]
    dfa = DFA(states, inputSize, final)

    dfa.addEdge(0, 1)
    dfa.addEdge(0, 2)
    dfa.addEdge(1, 3)
    dfa.addEdge(1, 2)
    dfa.addEdge(2, 2)
    dfa.addEdge(2, 4)
    dfa.addEdge(3, 1)
    dfa.addEdge(3, 4)
    dfa.addEdge(4, 3)
    dfa.addEdge(4, 2)


    # print "Transitions: "
    # print dfa.Transitions

    # logArray = dfa.computeDensityRecurrence(100050, 0, 50)
    # print logArray

    ComputeOnDfa(dfa)






def callWithDFAExample1():

    states = 6
    inputSize = 2
    final = [0, 0, 0, 1, 0, 1]
    dfa = DFA(states, inputSize, final)

    dfa.addEdge(0, 0)
    dfa.addEdge(0, 1)
    dfa.addEdge(1, 3)
    dfa.addEdge(1, 2)
    dfa.addEdge(2, 1)
    dfa.addEdge(2, 3)
    dfa.addEdge(3, 1)
    dfa.addEdge(3, 2)
    dfa.addEdge(4, 4)
    dfa.addEdge(4, 5)
    dfa.addEdge(5, 4)
    dfa.addEdge(5, 5)

    print "before Transitions: "
    print dfa.Transitions

    #states 5 and 6 are both useless, hopefully this function removes the two.
    dfa.PreProcessDFA()

    print "after Transitions: "
    print dfa.Transitions

    ComputeOnDfa(dfa)


def callWithDFAExample2():
    states = 8
    inputSize = 2
    final = [0, 0, 0, 0, 1, 0, 0, 1]
    dfa3 = DFA(states, inputSize, final)

    dfa3.addEdge(0, 1)
    dfa3.addEdge(0, 5)
    dfa3.addEdge(1, 2)
    dfa3.addEdge(1, 2)
    dfa3.addEdge(2, 3)
    dfa3.addEdge(2, 3)
    dfa3.addEdge(3, 4)
    dfa3.addEdge(3, 4)
    dfa3.addEdge(4, 1)
    dfa3.addEdge(4, 6)
    dfa3.addEdge(5, 6)
    dfa3.addEdge(5, 6)
    dfa3.addEdge(6, 7)
    dfa3.addEdge(6, 7)
    dfa3.addEdge(7, 5)
    dfa3.addEdge(7, 5)


    ComputeOnDfa(dfa3)


def ComputeOnDfa(dfa):

    task = int(input("1 for a graph of density trend; 2 for difference in time for single computation; 3 for more qualitative benchmarking."))
    if task == 1:
        wordLen = 100000
        lowLogValue = 99950
        highLogValue = 99999
        logArray = dfa.LogDensityRecurrence(wordLen, lowLogValue, highLogValue)
        print len(logArray)
        print logArray


    #matplotlib code
        index = []
        for i in range(highLogValue - lowLogValue +1):
            index.append(i + lowLogValue)

        plt.plot(index, logArray)

        plt.show()

    elif task == 2:
        wordLen = int(input("provide input length, below 1024 please:  "))
        #time both of these tasks, print out values
        matrixTimeStart = time.time()
        matrixDensity = dfa.ComputeDensityMatrix(wordLen +1)
        matrixTimeEnd = time.time()
        recurrenceTimeStart = time.time()
        recurrenceDensity = dfa.ComputeDensityRecurrence(wordLen)
        recurrenceTimeEnd = time.time()
        # print "Matrix Density: ", matrixDensity
        # print "Recurrence Density: ", recurrenceDensity
        print "\n"
        print "matrix time: " + str(matrixTimeEnd - matrixTimeStart)
        print "recurrence time: " + str(recurrenceTimeEnd - recurrenceTimeStart)
    elif task == 3:

        WordSet = []
        print "please provide inputs for matrix vs recurrence density computation."
        print "type 0 to finish entering data."
        doneEntering = True
        while(doneEntering):
            inputValue = int(input("provide input value: "))
            if inputValue == 0:
                doneEntering = False
            else:
                WordSet.append(inputValue)
        for i in WordSet:
            matrixTimeStart = time.time()
            matrixDensity = dfa.ComputeDensityMatrix(i +1)
            matrixTimeEnd = time.time()
            recurrenceTimeStart = time.time()
            recurrenceDensity = dfa.ComputeDensityRecurrence(i)
            recurrenceTimeEnd = time.time()
            print "Size: " + str(i) + " matrix: " + str(matrixTimeEnd - matrixTimeStart) + " recurrence: " + str(recurrenceTimeEnd - recurrenceTimeStart)

    else:
        exit(1)



def callWithNFAExample0():
    states = 3
    inputSize = 2
    Final = [0,0,1]
    nfa = NFA(states, inputSize, Final)

    nfa.addEdge(0,0,0)
    nfa.addEdge(0,0,1)
    nfa.addEdge(0,1,0)
    nfa.addEdge(1,2,0)
    nfa.addEdge(1,2,1)

    wordLen = 1000
    lowLogValue = 950
    highLogValue = 999
    logArray = nfa.LogDensityRecurrence(wordLen, lowLogValue, highLogValue)

    print logArray



def callWithNFAExample1():
    states = 12
    inputSize = 2
    Final = [0,0,0,0,1,0,0,0,0,0,1,0]
    nfa = NFA(states, inputSize, Final)


    nfa.addEdge(0,1,1)
    nfa.addEdge(0,1,0)
    nfa.addEdge(1,2,0)
    nfa.addEdge(1,2,1)
    nfa.addEdge(1,6,0)
    nfa.addEdge(2,3,0)
    nfa.addEdge(2,3,1)
    nfa.addEdge(2,11,0)
    nfa.addEdge(3,11,1)
    nfa.addEdge(3,4,0)
    nfa.addEdge(3,4,1)
    nfa.addEdge(4,5,0)
    nfa.addEdge(4,5,1)
    nfa.addEdge(5,2,0)
    nfa.addEdge(5,2,1)
    nfa.addEdge(6,7,0)
    nfa.addEdge(6,7,1)
    nfa.addEdge(7,8,0)
    nfa.addEdge(7,8,1)
    nfa.addEdge(8,9,0)
    nfa.addEdge(8,9,1)
    nfa.addEdge(9,10,0)
    nfa.addEdge(9,10,1)
    nfa.addEdge(10,11,1)
    nfa.addEdge(10,6,1)
    nfa.addEdge(10,6,0)
    nfa.addEdge(11,11,0)
    nfa.addEdge(11,11,1)
    # print nfa.Transitions


    wordLen = 100
    lowLogValue = 0
    highLogValue = 50
    logArray = nfa.LogDensityRecurrence(wordLen, lowLogValue, highLogValue)

    print logArray



def main():
    #states, transitions=None, inputSize, final:

    answer = int(input("0 for dfa or 1 for nfa: "))#1, 2, 3 or 0 for nfa
    if(answer == 0):
        answer = int(input("DFA 0, 1, or 2: "))
        if answer == 0:
            callWithDFAExample0()
        elif answer == 1:
            callWithDFAExample1()
        elif answer == 2:
            callWithDFAExample2()
        else:
            print " bad answer, exiting."
            exit(1)

    else:
        answer = int(input("NFA 0, or 1: "))
        if answer == 0:
            callWithNFAExample0()
        elif answer == 1:
            callWithNFAExample1()


    #   dfa.PreProcessDFA()

    return

main()
