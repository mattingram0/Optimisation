#    def slowPowerset(self, C, S):
#        if len(C) == 0:
#            connected = set()
#            for i in range(len(S)):
#                for j in range(i + 1, len(S)):
#                    #print('i: ' + str(i) + ', S[' + str(i) + '] = ' + str(list(S)[i]) + '. j: ' + str(j) + ', S[' + str(j) + '] = ' + str(list(S)[j]))
#                    if self.matrix[list(S)[i]][list(S)[j]] == 1:
#                        connected.add(list(S)[i]) 
#                        connected.add(list(S)[j])
#            if S == connected:
#                self.powerset.append(S)
#        else:
#            C1 = C.copy()
#            v = C1.pop()          
#            S1 = S.copy()
#            S2 = S.copy()
#            self.slowPowerset(C1, S1)
#            S2.add(v)
#            self.slowPowerset(C1, S2)
#
    #def subgraph(self, S, C, X):
    #   self.powerset.append(S)

   #    for v in C:
   #        if v in X: continue

   #        S1 = S.copy()
   #        S1.add(v)
   #        C1 = self.getNeighbours(v).difference(S)
   #        X1 = X.copy()
   #        self.subgraph(S1, C1, X1)
   #        X.add(v)
