import os, sys, itertools
from itertools import chain, combinations
from scipy import optimize as opt
from math import floor, inf
from pyomo.environ import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


class Graph:

    ''' Constructor for the Graph Class '''
    def __init__(self, name="", matrix=[], vertices=set(), edges=set(), cliques=[], powerset=[]):
        self.name = name            
        self.matrix = matrix        #NxN Adjacency Matrix
        self.vertices = vertices    #Set of all vertices
        self.edges = edges          #Set of all edges
        self.cliques = cliques      #List of all cliques
        self.powerset = powerset    #List of all subsets of graph
   
    '''
    Load a Graph From a File
    Args:
        filename - .txt file to load the graph from
    '''
    def loadGraph(self, filename):
        path = os.path.join(sys.path[0], filename) 
        with open(path, 'r') as f:
            readData = f.readline().split(' ')
            numVertices = int(readData[0])
            numEdges = int(readData[1])

            #0 the adjacency matrix
            matrix = [[0 for i in range(numVertices)] for j in range(numVertices)]

            #Insert a 1 at matrix[i][j] and matrix[j][i] if vertices i and j are linked by an edge 
            for i in range(numEdges):
                readData = f.readline().split(' ')
                i = int(readData[0])
                j = int(readData[1])
                matrix[i][j] = 1
                matrix[j][i] = 1

            self.matrix = matrix
        self.name = filename.split(".")[0]

    '''
    Find the Neighbours of a Vertex
    Args:
        v - a vertex represented as in integer
    Returns:
        A set containing all of the vertexes adjacent to v in the graph
    '''
    def getNeighbours(self, v):
        return set([i for i, x in enumerate(self.matrix[v]) if x == 1])

    '''
    Find all Cliques in the Graph (using a modified version of the Bron Kerbosch algorithm
    Args
        R - Set of vertices currently in the clique
        P - Set of candidate vertices
        X - Set of excluded vertices
    '''
    def bronKerbosch(self, R, P, X):
        if len(X) == 0:
            self.cliques.append(R)
        for v in P.copy():
            N = self.getNeighbours(v)
            R1 = set(R)
            R1.add(v)
            P1 = P.intersection(N)
            X1 = X.intersection(N)

            self.bronKerbosch(R1, P1, X1) 

            P.remove(v)

    '''
    Find all cliques in the graph. Effectively just a wrapper function for the recursive bronKerbosch function
    '''
    def findCliques(self):
        #Make a copy to ensure pass-by-value and not pass-by-reference
        self.bronKerbosch(set(), self.vertices.copy(), set())

    '''
    Find the set of all subsets of vertices of the graphs
    '''
    def findPowerset(self):
        self.powerset = list(chain.from_iterable(combinations(self.vertices, r) for r in range(len(self.vertices)+1)))

    '''
    Find all the subsets of a set of length between m and n
    Args
        s - the set to find the subsets of m
        m - the lower bound on the length of the subsets to find
        n - the upper bound on the length of the subsets to find
    Returns
        A list containing all the subsets of the specified lengths
    '''
    def getPowersetMN(self, s, m, n):
        return list(chain.from_iterable(combinations(s, r) for r in range(m, n)))

    '''
    Get a set of all the vertices in the graph
    '''
    def getVertices(self):
        self.vertices = set([i for i in range(len(self.matrix))])

    '''
    Get the objective function for the Fractional Clique Cover Number linear program
    Returns
        A list containing the coefficients of the variables in the objective function
    '''
    def getFCCNObjective(self):
        return [1 for i in range(len(self.cliques))]

    '''
    Get all of the constraints for the Fractional Clique Cover Number linear program
        Returns
            A list containing two lists:
                - list[0]: Contains the lists of all the coefficients of the LHS of the constraints
                - list[1]: Contains the RHS of every constraint
    '''
    def getFCCNConstraints(self):
        allLHS = []
        for v in self.vertices:
            lhs = [0 for i in range(len(self.cliques))]
            for c in self.cliques:
                if v in c:
                    lhs[self.cliques.index(c)] = 1
            allLHS.append(lhs)

        allRHS = [1 for i in range(len(allLHS))]
        return [allLHS, allRHS]

    '''
    Get the bounds of all the variables in the Fractional Clique Cover Number Linear Program
    '''
    def getFCCNBounds(self):
        return [(0, None) for i in range(len(self.cliques))]

    '''
    Create the linear programming model required to find the Fractional Clique Cover Number
    '''
    def createFCCNModel(self):
        objective = self.getFCCNObjective()
        constraints = self.getFCCNConstraints()
        bounds = self.getFCCNBounds()

        #Three inline functions used by Pyomo to create the model
        def eqConstraintRule(model, j):
            constraint = []
            return sum(model.x[k] for k, v in enumerate(constraints[0][j]) if v == 1) == 1

        def boundsRule(model, i):
            return bounds[i]
        
        def objectiveRule(model):
            return sum(model.x[i] for i in model.I)

        #Instantiate Model
        model = ConcreteModel()

        #Define Indexes
        model.I = RangeSet(0, len(self.cliques) - 1) #Variable Index
        model.J = RangeSet(0, len(self.vertices) - 1) #Constraint Index

        #Create Variables
        model.x = Var(model.I, domain=NonNegativeReals, bounds=boundsRule)

        #Create Objective Function
        model.obj = Objective(rule=objectiveRule, sense=minimize)

        #Create Constraints
        model.EqConstraint = Constraint(model.J, rule=eqConstraintRule)
        
        #Write the model to a file
        model.write(self.name + '_FCCN.mps', format='mps')

    '''
    Get the objective function for the Shannon Entropy linear program
    Returns
        A list containing the coefficients of the variables in the objective function
    '''
    def getSEObjective(self):
        #List of 0s except for the final elementr
        return [-1 * floor(i/(len(self.powerset) - 1)) for i in range(len(self.powerset))]

    '''
    Get the constraints for the Shannon Entropy linear program
    Returns
        A list containing four lists:
            - list[0]: Contains the lists of all the coefficients of the LHS of the equality constraints
            - list[1]: Contains the RHS of all the equality constraints
            - list[2]: Contains the lists of all the coefficients of the LHS of the less than constraints
            - list[3]: Contains the RHS of all the less than constraints
            
    '''
    def getSEConstraints(self):
        eqlLHS, eqlRHS, ltLHS, ltRHS = [], [], [], []

        #Empty set equality
        eqlLHS.append([floor(1/i) for i in range(1, len(self.powerset) + 1)])
        eqlRHS.append(0)

        #Monotone increasing and neighbourhood inequalities
        for t in self.powerset:
            if len(t) <= 1: continue
            subpowerset = self.getPowersetMN(t, len(t) - 1, len(t))
            for s in subpowerset:
                constraint = [0 for i in range(len(self.powerset))]
                #NOTE THAT ALL GREATER THAN CONSTRAINTS HAVE BEEN INVERTED TO BE LESS THAN CONSTRAINTS
                constraint[self.powerset.index(t)] = -1
                constraint[self.powerset.index(s)] = 1
               
                difference = set(t).difference(set(s))
                dif = difference.copy()
                vertex = difference.pop()
                if set(t) == dif.union(self.getNeighbours(vertex)):
                    eqlLHS.append(constraint)
                    eqlRHS.append(0)
                else:
                    ltLHS.append(constraint)
                    ltRHS.append(0)

        #Union and intersection constraint
        for c in combinations(self.powerset, 2):
            constraint = [0 for i in range(len(self.powerset))]
            pair = list(c)
            s = set(pair[0])
            t = set(pair[1])
            union = s.union(t)
            intersection = s.intersection(t)

            if s.issubset(t) or t.issubset(s): continue 

            #NOTE THAT ALL GREATER THAN CONSTRAINTS HAVE BEEN INVERTED TO BE LESS THAN CONSTRAINTS
            constraint[self.powerset.index(tuple(s))] = -1
            constraint[self.powerset.index(tuple(t))] = -1
            constraint[self.powerset.index(tuple(union))] = 1
            constraint[self.powerset.index(tuple(intersection))] = 1
            ltLHS.append(constraint)
            ltRHS.append(0)

        return [eqlLHS, ltLHS, eqlRHS, ltRHS]

    '''
    Get the bounds of all the variables in the Shannon Entropy Linear Program
    '''
    def getSEBounds(self):
        #All variables have to be greater than 0
        bounds = [(0, None) for i in range(len(self.powerset))]
        
        for i in range(len(self.vertices)):
            #Variables corresponding to single-vertex subsets are upper-bounded by 1
            bounds[i + 1] = (0, 1)
        return bounds
    
    '''
    Debugging function used to print out the formulated Linear Program
    '''
    def printLP(self, constraints, objective, bounds):
        print("Objective Function:")

        obj = ''
        for i in range(len(self.powerset)):
            if objective[i] == 0: 
                continue
            elif objective[i] == 1:
                obj = obj + ' + ' + 'x' + str(i) 
            else:
                obj = obj + ' - ' + 'x' + str(i) 

        print(obj) 
        print()
        print("Variable Allocations:")
        counter = 0
        for s in self.powerset:
            print('x' + str(counter) + ':')
            print(s)
            counter+=1

        print()
        print("Constraints:")

        eqlLHS = constraints[0]
        eqlRHS = constraints[2]
        upLHS = constraints[1]
        upRHS = constraints[3]

        for e in eqlLHS:
            constraint = ''
            for i in range(len(self.powerset)):
                if e[i] == 0: 
                    continue
                elif e[i] == 1:
                    constraint = constraint + ' + ' + 'x' + str(i) 
                else:
                    constraint = constraint + ' - ' + 'x' + str(i) 
                    
            constraint += ' = ' + str(eqlRHS[eqlLHS.index(e)])
            print(constraint)

        for u in upLHS:
            constraint = ''
            for i in range(len(self.powerset)):
                if u[i] == 0: 
                    continue
                elif u[i] == 1:
                    constraint = constraint + ' + ' + 'x' + str(i) 
                else:
                    constraint = constraint + ' - ' + 'x' + str(i) 
            constraint += ' <= ' + str(upRHS[upLHS.index(u)])
            print(constraint)

        print()
        print('Bounds:')

        for i in range(len(self.powerset)):
            print('x' + str(i) + ' is in [' + str(bounds[i][0]) + ', ' + str(bounds[i][1] if bounds[i][1] != None else inf) + ']')

        print()

    '''
    Create the linear programming model required to find the Shannon Entropy
    '''
    def createSEModel(self):
        objective = self.getSEObjective()
        constraints = self.getSEConstraints()
        bounds = self.getSEBounds()

        def eqConstraintRule(model, j):
            constraint = []
            for idx, val in enumerate(constraints[0][j]):
                if val == 1:
                    constraint.append(idx)
                elif val == -1:
                    constraint.append(-1 * idx)

            return sum(model.x[0] if j == 0 else int(j/abs(j)) * model.x[abs(j)] for j in constraint) == 0

        def ltConstraintRule(model, k):
            constraint = []
            for idx, val in enumerate(constraints[1][k]):
                if val == 1:
                    constraint.append(idx)
                elif val == -1:
                    constraint.append(-1 * idx)

            return sum(model.x[0] if j == 0 else int(j/abs(j)) * model.x[abs(j)] for j in constraint) <= 0

        def boundsRule(model, i):
            return bounds[i]

        #Instantiate Model
        model = ConcreteModel()

        #Define Indexes
        model.I = RangeSet(0, len(self.powerset) - 1) #Variable Index
        model.J = RangeSet(0, len(constraints[0]) - 1) #Equality Constraint Index
        model.K = RangeSet(0, len(constraints[1]) - 1) #Less Than Constraint Index

        #Create Variables
        model.x = Var(model.I, domain=NonNegativeReals, bounds=boundsRule)

        #Create Objective Function
        model.obj = Objective(expr = (-1) * model.x[len(self.powerset) - 1])

        #Create Constraints
        model.EqConstraint = Constraint(model.J, rule=eqConstraintRule)
        model.LtConstraint = Constraint(model.K, rule=ltConstraintRule)

        #Write Model to mps file 
        model.write(self.name + '_SE.mps', format='mps')
    
def main():
    graph = Graph()
    graph.loadGraph("sys.argv[1]")
    graph.getVertices()
    graph.findCliques()
    graph.findPowerset()

    graph.createFCCNModel()
    graph.createSEModel()

if __name__ == "__main__":
    main()
