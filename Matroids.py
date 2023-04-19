#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools
import sympy as sym
from rustworkx.visualization import graphviz_draw
import rustworkx as rx

class AbstractMatroid:
    def __init__(self,baseSet,independentSets):
        self.E = baseSet
        self.I = independentSets
        
    @property
    def cardinality(self):
        return len(self.E)

    @property
    def totalrank(self):
        return max([len(t) for t in self.I])
    
    @property
    def bases(self):
        n = self.totalrank
        subs = self.allSubsets
        result = set()
        for sub in subs:
            if self.isBasis(sub):
                result.add(frozenset(sub))
        return result
    
    @property
    def flats(self):
        result = set()
        for sub in self.I:
            result.add((self.rank(sub),frozenset(self.closure(sub))))
        return result
    
    @property
    def allSubsets(self):
        return self.subsets(self.cardinality)
    
    @property
    def characteristicPolynomial(self):
        t = sym.Symbol('t')
        result = 0
        for sub in self.allSubsets:
            result += (-1)**len(sub)*t**(self.totalrank - self.rank(sub))
        return result
    
    @property
    def reducedCharacteristicPolynomial(self):
        t = sym.Symbol('t')
        return sym.simplify(self.characteristicPolynomial/(t-1))
    
    @property
    def beta(self):
        t = sym.Symbol('t')
        return (-1)**(self.totalrank-1)*self.reducedCharacteristicPolynomial.subs(t,1)
    
    @property
    def lattice(self):
        graph = rx.PyDiGraph(check_cycle=True)
        graph.add_nodes_from(list(self.flats))
        nodeIndices = graph.node_indices()
        for i in nodeIndices:
            dim = self.totalrank
            n1 = graph[i][0]
            flat1 = graph[i][1]
            if n1 < dim:
                for j in nodeIndices:
                    n2 = graph[j][0]
                    if n2 > n1:
                        flat2 = graph[j][1]
                        if flat1.issubset(flat2):
                            graph.add_edge(i,j,n2-n1)
        return graph
    
    @property
    def reducedLattice(self):
        graph = rx.PyDiGraph(check_cycle=True)
        graph.add_nodes_from(list(self.flats))
        nodeIndices = graph.node_indices()
        for i in nodeIndices:
            dim = self.totalrank
            n1 = graph[i][0]
            flat1 = graph[i][1]
            if n1 < dim:
                for j in nodeIndices:
                    n2 = graph[j][0]
                    if n2 - n1 == 1:
                        flat2 = graph[j][1]
                        if flat1.issubset(flat2):
                            graph.add_edge(i,j,None)
        return graph
    
    @property
    def terminalFlats(self):
        graph = self.lattice
        result = [0,0]
        for i in graph.node_indices():
            if graph[i][0] == 0:
                result[0] = i
            elif graph[i][0] == self.totalrank:
                result[1] = i
        return result
    
    @property
    def chains(self):
        graph = self.lattice
        terminal = self.terminalFlats
        lchains = rx.digraph_all_simple_paths(graph,terminal[0],terminal[1])
        return [[graph[i] for i in chain] for chain in lchains]
    
    @property
    def chainIndices(self):
        graph = self.lattice
        terminal = self.terminalFlats
        return rx.digraph_all_simple_paths(graph,terminal[0],terminal[1])
    
    def restrict(self,s):
        self.E = s
        self.I = set([t.intersection(self.E) for t in self.I])
        
    def dualize(self):
        ind = set()
        n = self.totalrank
        for sub in self.allSubsets:
            if self.rank(sub) == n:
                ind.add(frozenset(self.E.difference(sub)))
        self.I = ind
        
    def dual(self):
        temp = AbstractMatroid(self.E,self.I)
        return temp.dualize()
    
    def contract(self,s):
        if s != set():
            ind1 = self.I
            temp = AbstractMatroid(self.E,ind1)
            temp.restrict(s)
            self.restrict(self.E.difference(s))
            bases = temp.bases
            ind2 = self.allSubsets
            result = set()
            for sub in ind2:
                for basis in bases:
                    if {sub.union(basis)}.intersection(ind1) != set():
                        result = result.union(set([frozenset(sub)]))
            self.I = result
        
    
    def rank(self,s):
        temp = AbstractMatroid(self.E,self.I)
        temp.restrict(s)
        return temp.totalrank
    
    def subsets(self,r):
        sub = {frozenset({})}
        for i in range(1,r+1):
            sub = sub.union(set([frozenset(i) for i in itertools.combinations(self.E,i)]))
        return sub
    
    def isBasis(self,s):
        if self.isIndependent(s):
            if self.rank(s) == self.totalrank:
                return True
            else:
                return False
        else:
            return False
        
    def isIndependent(self,s):
        if s in self.I:
            return True
        else:
            return False
    
    def closure(self,s):
        n = self.rank(s)
        result = set(s)
        complement = self.E
        complement = complement.difference(result)
        for i in complement:
            if self.rank(s.union({i})) == n:
                result.add(i)
        return result
    
    def drawLattice(self):
        return graphviz_draw(self.reducedLattice, node_attr_fn = lambda node:{"label":str(set(node[1]))})
    
    def cones(self,k):
        lchains = self.chains
        result = []
        for chain in lchains:
            if len(chain) == k+2:
                result.append(chain)
        return result
    
    def csmWeight(self,chain):
        result = (-1)**(self.totalrank - len(chain) + 1)
        for i in range(0,len(chain)-1):
            temp = AbstractMatroid(self.E,self.I)
            temp.restrict(chain[i+1][1])
            temp.contract(chain[i][1])
            result = result*temp.beta
        return result
            
    
    def relabel(self,f):
        self.E = set([f(i) for i in self.E])
        self.I = set(frozenset([f(i) for i in sub]) for sub in self.I)
        
    def chowRing(self):
        chowVars = dict()
        flats = set([flat[1] for flat in self.flats if not (flat[1].intersection(self.E) == self.E or flat[1] == set({}))])
        for flat in flats:
            chowVars[flat] = sym.Symbol('x_'+str(set(flat)))
        baseRing = sym.QQ.old_poly_ring(*list(chowVars.values()))
        pairsIdeal = []
        for flat1 in flats:
            for flat2 in flats.difference({flat1}):
                if not (flat1.issubset(flat2) or flat2.issubset(flat1)):
                    pairsIdeal.append(chowVars[flat1]*chowVars[flat2])
        groundElements = list(self.E)
        differenceIdeal = []
        while len(groundElements) > 1:
            i = groundElements.pop(0)
            for j in groundElements:
                contains_i = [chowVars[flat] for flat in flats if flat.intersection({i}) != set({})]
                contains_j = [chowVars[flat] for flat in flats if flat.intersection({j}) != set({})]
                differenceIdeal.append(sum(contains_i) - sum(contains_j))
        idealSum = (pairsIdeal + differenceIdeal)
        return baseRing.quotient_ring(idealSum)
    
    def __add__(self,o):
        if self.E.intersection(o.E) == set([]):
            base = self.E.union(o.E)
            ind = set([frozenset(sub1.union(sub2)) for sub1,sub2 in itertools.product(self.I,o.I)])
        else:
            base = set(itertools.product(self.E,set([0]))).union(set(itertools.product(o.E,set([1]))))
            ind1 = set([ frozenset([(i,0) for i in sub]) for sub in self.I])
            ind2 = set([ frozenset([(i,1) for i in sub]) for sub in o.I])
            ind = set([frozenset(sub1.union(sub2)) for sub1,sub2 in itertools.product(ind1,ind2)])
        return AbstractMatroid(base,ind)
    
    def __str__(self):
        result = f"Base set: \n{self.E}, \nIndependent sets:"
        for sub in self.I:
            result+="\n"+str(sub)
        return result

def uniformMatroid(n,r):
    base = frozenset(range(0,n))
    ind = {frozenset({})}
    for i in range(1,r+1):
        ind = ind.union(set([frozenset(i) for i in itertools.combinations(base,i)]))
    return AbstractMatroid(base,ind)

def starMatroid(n):
    m = uniformMatroid(n+1,2)
    for i in range(1,n):
        for j in range(i,n+1):
            m.I.add(frozenset({0,i,j}))
    return m

