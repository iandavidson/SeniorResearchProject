import numpy as np
import sys
from collections import defaultdict
import copy as copy
import matplotlib.pyplot as plt
import time as time




###############################################################################
# NFA class: to represent a regular language for usage for the data model : Unambiguous NFA #
#                                                     #
# members:
#   States -> number of states in instance
#   Transitions -> [States -> [InputSize -> [set of transitions]]] 3d list, list can be passed in or addEdge can be used to add edges in graph.(__init__)
#   Final -> list size States, that have either 0 or 1, whether or not they are accepting.
#   InputSize -> the size of the input language. if we are using binary inputs, the value would be 2.
# 
#  Mutator methods:
#     addEdge -> used if we start with a blank transition table, we can fill in items one at a time
#     PreProcessNFA -> used to fix issues with the given instance, usually should be called after creation.
#       FixNullDeltas -> will add self loops for transitions values that are null (that way we mitigate edge cases when transitions are null).
#       InitialBreadthFirstSearch -> method used to go through one iteration of a BFS given starting state. Visited Vector is returned, conveying which states are reachable from the starting state.
# 
# 
#  Operations of Class:
#     delta -> takes in a state and input language index, and returns the state index
#     LogDensityRecurrence -> returns a list of density computations given arguements.
#     ComputeDensityRecurrence -> returns a float for input length and DFA instance, using a Theta(n) recurrence relation method.
#     ComputeDensityMatrix -> returns a float for input length and DFA instance, using a Matrix Method method.
#       divideByConstantHelper -> computes matrix multiplication in Theta(logn)
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

    #function used to add an edge directed from stateIdx: <fromState> -> <endState>
    def addEdge(self, fromState, endState, edgeValue):
        self.Transitions[fromState][edgeValue].append(endState)



    # We want to see if there are any unreachable states from the start
    # and remove them, we check at multiple
    def PreProcessNFA(self):

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
                #input0 /= len(delta0)
                input1 = 0.0
                delta1 = self.delta(j, 1)
                for l in range(len(delta1)):
                    input1 = float(input1 + old[delta1[l]])
                #input1 /= len(delta1)

                #if input0 or input1 is  greater than 1 at any point, then we can say the nfa is NOT unambiguous


                #current[j] = float(input0 + input1)/2.0 #force a cast to float type
                # if len(delta0)+ len(delta1) > 1:
                #     current[j] = float(input0 + input1)/(len(delta0) + len(delta1)) #force a cast to float type
                # else:
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
        #this function will only compute density at 
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


