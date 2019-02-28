#code ripped from here: https://www.geeksforgeeks.org/strongly-connected-components/

# Python implementation of Kosaraju's algorithm to print all SCCs
from collections import defaultdict



Transitions = [[2,1],[3,2],[4,2],[1,4],[3,2]]


#This class represents a directed graph using adjacency list representation
class Graph:

	def __init__(self,vertices):
		self.V= vertices #No. of vertices
		self.graph = defaultdict(list) # default dictionary to store graph
		#member for strongly connected components


	# function to add an edge to graph
	def addEdge(self,u,v):
		self.graph[u].append(v)

	# A function used by DFS
	def DFSUtil(self,v,visited):
		# Mark the current node as visited and print it
		visited[v]= True
		print v,
		#Recur for all the vertices adjacent to this vertex
		for i in self.graph[v]:
			if visited[i]==False:
				self.DFSUtil(i,visited)


	def fillOrder(self,v,visited, stack):
		# Mark the current node as visited
		visited[v]= True
		#Recur for all the vertices adjacent to this vertex
		for i in self.graph[v]:
			if visited[i]==False:
				self.fillOrder(i, visited, stack)
		stack = stack.append(v)


	# Function that returns reverse (or transpose) of this graph
	def getTranspose(self):
		g = Graph(self.V)

		# Recur for all the vertices adjacent to this vertex
		for i in self.graph:
			for j in self.graph[i]:
				g.addEdge(j,i)
		return g



	# The main function that finds and prints all strongly
	# connected components
	def printSCCs(self):

		stack = []
		# Mark all the vertices as not visited (For first DFS)
		visited =[False]*(self.V)
		# Fill vertices in stack according to their finishing
		# times
		for i in range(self.V):
			if visited[i]==False:
				self.fillOrder(i, visited, stack)

		# Create a reversed graph
		gr = self.getTranspose()

		# Mark all the vertices as not visited (For second DFS)
		visited =[False]*(self.V)

		# Now process all vertices in order defined by Stack
		while stack:
			i = stack.pop()
			if visited[i]==False:
				gr.DFSUtil(i, visited)
				print""




def main():
        # lets use this set of transitions: global Transitions
        # Create a graph given in the above diagram
        g = Graph(len(Transitions))

		# create an instance of a graph and add in edges from transition tableself.
		#note the graph class stores these in the form of an adjacency list.
        for i in range(len(Transitions)):
                for j in range(len(Transitions[i])):
                        g.addEdge(i, Transitions[i][j])

        print ("Following are strongly connected components " +
               "in given graph")
        g.printSCCs()
		#we also need the components stored in some way

        #  Compute Period: (A, s) =>
        #     A => an n by n transition table holding probabilities between 0 and 1
        #     s => intial spread -> START[] -> can be literally taken.
        #       find subchains A[1] -> A[p]  induced by terminal components with periodicy of l[1] -> l[p], we say there are p subchains of A
        #      {might be wrong part}
        #     For i=0 -> p; j=0 -> "l[i]-1":
        #         do => alpha[i][j] = SumofAcceptingStates(Compute-Limit(A^(l[i]), A^(j)s))
        #     return alpha //a set of |Q| lists
        #
        #  Compute Limit: (B, s) =>
        #    B => a transition matrix
        #     Compute Terminal
        #
        #
        #
        #
        #






main()
