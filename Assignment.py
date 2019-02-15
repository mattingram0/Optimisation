import os, sys, itertools
from itertools import chain, combinations
from scipy import optimize as opt
from math import floor, inf
from pyomo.environ import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


class Graph:

    def __init__(self, name="", matrix=[], vertices=set(), edges=set(), cliques=[], powerset=[]):
        self.name = name
        self.matrix = matrix        #NxN Adjacency Matrix
        self.vertices = vertices    #Set of all vertices
        self.edges = edges          #Set of all edges
        self.cliques = cliques      #List of all cliques
        self.powerset = powerset    #List of all subsets of graph
    
    def loadGraph(self, filename):
        path = os.path.join(sys.path[0], filename) 
        with open(path, 'r') as f:
            readData = f.readline().split(' ')
            numVertices = int(readData[0])
            numEdges = int(readData[1])

            #0 the adjacency matrix
            matrix = [[0 for i in range(numVertices)] for j in range(numVertices)]

            #Insert a 1 at matrix[i][j] and matrix[j][i] is vertices i and j are linked by an edge 
            for i in range(numEdges):
                readData = f.readline().split(' ')
                i = int(readData[0])
                j = int(readData[1])
                matrix[i][j] = 1
                matrix[j][i] = 1

            self.matrix = matrix
        self.name = filename.split(".")[0]

    def getNeighbours(self, v):
        return set([i for i, x in enumerate(self.matrix[v]) if x == 1])

    #Algorithm to find all the maximal cliques in the graph
    #Modified version of the bronKerbosch algorithm, which finds all cliques
    #and not only the maximal ones
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

    def findCliques(self):
        #Make a copy to ensure pass-by-value and not pass-by-reference
        self.bronKerbosch(set(), self.vertices.copy(), set())

    def getVertices(self):
        self.vertices = set([i for i in range(len(self.matrix))])

    '''Returns a zero-list with length equal to the number of cliques'''
    def getFCCNObjective(self):
        return [1 for i in range(len(self.cliques))]

    '''Returns a list containing two list:
        - The first contains all the lists of the coefficients of the LHS of the constraints
        - The second contains all the RHS of the constraints'''
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

    '''Returns the Fractional Clique Cover Number of the Graph, using Linear Programming'''
    def findFCCN(self):
        constraints = self.getFCCNConstraints()
        objective = self.getFCCNObjective()
        result = opt.linprog(objective, A_eq=constraints[0], b_eq=constraints[1])
        print("----- Fractional Clique Cover Number -----")
        print("Value: " + str(result.fun))
        output = "Found At: "
        for i in range(len(self.cliques) - 1):
            output += "x[" + str(i) + "] = " + str(result.x[i]) + ", "
        output += "x[" + str(len(self.cliques) - 1) + "] = " + str(result.x[len(self.cliques) - 1]) 
        print(output)

    def getSEObjective(self):
        #List of 0s except for the final elementr
        return [-1 * floor(i/(len(self.powerset) - 1)) for i in range(len(self.powerset))]

    def findPowerset(self):
        self.powerset = list(chain.from_iterable(combinations(self.vertices, r) for r in range(len(self.vertices)+1)))

    '''Gets all the subsets of a set (s) of length between m and n'''
    def getPowersetMN(self, s, m, n):
        return list(chain.from_iterable(combinations(s, r) for r in range(m, n)))

    def getSEConstraints(self):
        #This may be wrong - TEST
        eqlLHS, eqlRHS, ltLHS, ltRHS = [], [], [], []

        #Empty set equality
        eqlLHS.append([floor(1/i) for i in range(1, len(self.powerset) + 1)])
        eqlRHS.append(0)

        #Single vertex inequalities
        #ltLHS.extend([0**abs(i) for i in range(0 - j, len(self.powerset) + 1 - j)] for j in range(0, len(self.vertices)))
        #ltRHS.extend([1 for i in range(0, len(self.vertices))])

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
                #print(len(difference))
                dif = difference.copy()
                vertex = difference.pop()
                #difference.pop() will return the only element in the list - needs to be an int to pass to getNeighbours
                if set(t) == dif.union(self.getNeighbours(vertex)):
                    #print('T:')
                    #print(set(t))
                    #print('S:')
                    #print(set(s))
                    #print('set(t).union(self.getNeighbours(vertex)):')
                    #print(dif.union(self.getNeighbours(vertex)))
                    #print('Difference:')
                    #print(dif)
                    #print('Neighbours:')
                    #print(self.getNeighbours(vertex)) 
                    #print()
                    eqlLHS.append(constraint)
                    eqlRHS.append(0)
                else:
                    ltLHS.append(constraint)
                    ltRHS.append(0)

        #print('No. Equality Constraints')
        #print(len(eqlLHS))
        #print('No. Upper Bound Constraints')
        #print(len(ltLHS))

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

        #print('No. Upper Bound Constraints')
        #print(len(ltLHS))

        return [eqlLHS, ltLHS, eqlRHS, ltRHS]


    def getSEBounds(self):
        #All variables have to be greater than 0
        bounds = [(0, None) for i in range(len(self.powerset))]
        
        for i in range(len(self.vertices)):
            #Variables corresponding to single-vertex subsets are upper-bounded by 1
            bounds[i + 1] = (0, 1)
        return bounds

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

    def plotGraph(self):
        G = nx.from_numpy_matrix(np.array(self.matrix))
        nx.draw(G, with_labels=True)
        plt.show()

    def convertConstraints(self, constraints):
        eqlCnst = constraints[0]
        ltCnst = constraints[1]
        
        eqlCnstCon = []
        ltCnstCon = []

        for e in eqlCnst:
            eqlCnstCon.append([i for i in range(len(e)) if e[i] == 1])

        for l in ltCnst:
            ltCnstCon.append([i for i in range(len(l)) if l[i] == 1])

        return [eqlCnstCon, ltCnstCon]
        


    def findSEPyomo(self):
        objective = self.getSEObjective()
        constraints = self.getSEConstraints()
        bounds = self.getSEBounds()

        #Inline functions are nasty, but super useful here
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

        model.write(self.name + '.mps', format='mps')

        #Pretty Print Model
        #solver = SolverFactory('glpk')
        #solver.options['exact'] = 1
        #results = solver.solve(model)
        #results = SolverFactory('glpk').solve(model, options="exact=1")
        ##print('Write')
        #results.write()
#       # print('Pretty Print')
#        model.pprint()
#        print('Display')
#        model.x.display()
#
    def findSE(self):
        constraints = self.getSEConstraints()
        objective = self.getSEObjective()
        bounds = self.getSEBounds()
        #self.printLP(constraints, objective, bounds)
        #self.plotGraph() 
       # print('Equality Constraints (1 + 4 = 5)')
       # print(constraints[0])
       # print(len(constraints[0]))
       # print()
       # print('Upper Bound Constraints (28 - 4 + 120 - 16 - 61')
       # print(constraints[1])
       # print(len(constraints[1]))
       # #print(constraints[1])
        #print(constraints[2])
        #print(constraints[3])
        #print()
        #print('Number Of Constraints:')
        print(len(constraints[0]) + len(constraints[1]))
        #print(opt.show_options(solver="linprog"))
        #result = opt.linprog(objective, A_eq=constraints[0], b_eq=constraints[2], A_ub=constraints[1], b_ub=constraints[3], bounds=bounds, method='simplex', options={'tol':1e-8})
        #result = opt.linprog(objective, A_eq=constraints[0], b_eq=constraints[2], A_ub=constraints[1], b_ub=constraints[3], bounds=bounds, method='interior-point', options={'tol':1e-13})

        print("----- Shannon Entropy -----")
        print("Value: " + str(abs(result.fun)))
#       output = "Found At: "
#       for i in range(len(self.powerset) - 1):
#            output += "x[" + str(i) + "] = " + str(result.x[i]) + ", "
#        output += "x[" + str(len(self.powerset) - 1) + "] = " + str(result.x[len(self.powerset) - 1]) 
#        print(output)
#
    
def main():
    graph = Graph()
    graph.loadGraph("complete3.txt")
    graph.getVertices()
    graph.findCliques()
    #graph.findFCCN()
    graph.findPowerset()
    #print(len(graph.powerset))
    #graph.findSE()
    graph.findSEPyomo()

if __name__ == "__main__":
    main()
