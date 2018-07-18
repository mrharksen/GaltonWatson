import numpy as np
import random
from scipy.special import zeta
import matplotlib.pyplot as plt
import networkx as nx
import pygraphviz as pgv
from networkx import Graph


'''
	Input:
			n           integer    large integer.
			p0          float      in (0,1) corresponding to probability of having 0 children.
			alpha       float      in (1,2) corresponding to exponent of distribution.
			gamma       float      in (0,1) and greater than alpha - 1.
			delta       float      small floating number.
			samplesize  integer    an integer corresponding to how many times
			draw        boolean    if True then the Galton Watson tree is drawn.

    Description:
	        We have two conditions:

			      p0 + p1 + p2 + ... = 1
				  0*p0 + 1*p1 + ... = 1

		where p0 is a given probability corresponding to a leaf not
		having any children, p1 is found by the conditions above.
        The rest of the probabilities follow pi ~ c*i^{-(alpha+1)} where alpha
		is the given input above. The constant c is found by the conditions above.

        We start by generating a Galton Watson tree from the distribution which the
		above values determine. From these values we can calculate the
'''
def GaltonWatson(n = 1000000000, p0 = 0.5, alpha = 1.5, gamma = 0.75, delta = 0.0001, samplesize = 5000, draw = True):
	c = findConstant(p0,alpha)
	p1 = findp1(p0,alpha,c)
	dist = Distribution(p0, p1, alpha, c)
	sumList = [0 for i in xrange(samplesize)]
	for i in xrange(samplesize):
		G = GaltonWatsonGenerator(dist)
		sumList[i] = sumLeaf(G, n, alpha, gamma, delta)
		if draw:
			H = nx.nx_agraph.to_agraph(G)
			H.layout(prog='dot')
			H.draw('test' + str(i) + '.png')
	return sumList
'''
	Creates the Galton Watson tree from the given distribution dist
'''
def	GaltonWatsonGenerator(dist):
	G = nx.Graph()
	G.add_node(0)
	nodeCount = 1
	bool = True
	covered = [False for i in xrange(10000)]
	while bool:
		bool = False
		if G.number_of_nodes() == 1:
			numChildren = getNumberOfChildren(dist)
			if (numChildren > 0):
				bool = True
			covered[0] = True
			while(numChildren > 0):
				G.add_node(nodeCount)
				G.add_edge(0,nodeCount)
				nodeCount += 1
				numChildren -= 1
		else:
			leaves = [node for node in G if G.degree(node) == 1]
			if nodeCount > len(covered):
				tmp = [False for i in xrange(2*nodeCount)]
				for i in xrange(len(covered)):
					tmp[i] = covered[i]
				del covered
				covered = tmp
			for w in leaves:
				if(covered[w] == True):
					pass
				else:
					covered[w] = True
					numChildren = getNumberOfChildren(dist)
					if (numChildren > 0):
						bool = True
					while (numChildren > 0):
						G.add_node(nodeCount)
						G.add_edge(w,nodeCount)
						nodeCount += 1
						numChildren -= 1
	return G

'''
	Calculates the function

	  sup_u sum_{v <= u} (d^+(v)/n^{1/alpha})^gamma 1_{d^+(v) <= delta*n^{1/alpha}}
'''
def sumLeaf(G, n, alpha, gamma, delta):
	mysum = 0
	leaves = [node for node in G if G.degree(node) == 1]
	for u in leaves:
		tempsum = 0
		pathToRoot = nx.shortest_path(G,source=u,target=0)
		for v in pathToRoot:
			if (G.degree(v)-1 <= delta*np.power(n,1/alpha)):
				tempsum += np.power((G.degree(v)-1)/np.power(n,1/alpha),gamma)
		if mysum < tempsum:
			mysum = tempsum
	return mysum

'''
    Finds the constant c by the conditions above.
'''
def findConstant(p0, alpha):
	return p0/(zeta(alpha)-zeta(alpha+1))

'''
    Finds the probability p1 by the conditions above.
'''
def findp1(p0, alpha, c):
	return 1+c-c*zeta(alpha)


'''
    Generates the distribution function of the corresponding probabilities
	p0, p1, ... such that sum pi > 0.99999999 (which can be altered).
'''
def Distribution(p0, p1, alpha, c):
	invDist = []
	invDist.append(p0)
	invDist.append(p0+p1)
	cumsum = p0+p1
	j = 1
	while(0.99999999 > cumsum):
		cumsum += c*np.power((j+1),-(alpha+1))
		invDist.append(cumsum)
		j += 1
	return invDist

'''
    Finds the size of the step which the random walk takes next.
'''
def getNumberOfChildren(dist):
	r = np.random.uniform()
	for x in dist:
		if r < x:
			return dist.index(x)
	return len(dist)-1

def main():
	n = 10000000000
	p0 = 0.2
	alpha = 1.5
	gamma = 0.75
	delta = 0.000001
	samplesize = 100
	sumList = GaltonWatson(n, p0, alpha, gamma, delta, samplesize, draw = False)
	for x in sumList:
		print x

if __name__ == "__main__":
	main()
