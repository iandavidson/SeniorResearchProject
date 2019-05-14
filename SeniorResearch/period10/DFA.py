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

        #if we have any transitions that are null, we need to add state then direct all edges to it.
        self.FixNullDeltas()

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


    def FixNullDeltas(self):
        #mutator function, helper for PreProcessDFA()
        #we need to look thorugh all of the transitions to see if there are 2 coming out of each state,
            #if we dont find any 'null' transitions, self object isnot altered.
            #else we add in another state that can be considered fail, this state will have all originally null deltas as incoming edges,
        needFailState = False
        for i in range(self.States):
            if len(self.Transitions[i]) < self.InputSize:
                needFailState = True
                dif = self.InputSize - len(self.Transitions[i])
                for j in range(dif):
                    self.addEdge(i, self.States)
            # for j in range(self.InputSize):
            #     if self.Transitions[i][j] ==

        if needFailState:
            deadState = [self.States, self.States]
            self.Transitions.push(deadState)
            self.States += 1
            self.FINAL = self.FINAL + [0]





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



    def getTranspose(self):
        #when a dfa is reversed, it becomes an unambiguous nfa
        #that means we will be returning an instance of NFA class.
        nfa = NFA(self.States, self.InputSize, self.Final)

        for i in range(States):
            for j in range(len(self.Transitions[i])):
                nfa.addEdge(self.Transitions[i][j], i)

        return nfa


    def SCCSearchHelper(self, index, visited, stack):
        visited[i] = True
        for i in range(len(self.Transitions[index])):
            if not visited[i]:
                self.SCCSearchHelper(index, visited, stack)
        stack = stack.append(v)

    def findSCCs(self):
        #unimplemented

        #make stack
        stack = []

        visited = [False] * self.States

        for i in range(len(self.Transitions)):
            if visited[i] == False:
                self.SCCSearchHelper(i, visted, stack)

        transpose = self.getTranspose()

        revisited = [False] * self.States

        componentIndex = 0
        while stack:
            index = stack.pop()
            if revisited[index] == False:
                self.StronglyCCs.append([])
                transpose.SCCSearchPrint(index, revisited, componentIndex )
                componentIndex += 1

        for j in range(len(transpose.StronglyCCs)):
            self.StronglyCCs.append(copy.deepcopy(transpose.StronglyCCs[j]))

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
