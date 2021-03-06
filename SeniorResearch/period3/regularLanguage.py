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


def callWithExample1():
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


    print "Transitions: "
    print dfa.Transitions


def callWithExample2():

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

def main():
    #states, transitions=None, inputSize, final:

    answer = int(input("call with dfa 1, or 2: "))
    if(answer == 1):
        callWithExample1()
    elif(answer == 2):
        callWithExample2()


    #   dfa.PreProcessDFA()

    return

main()
