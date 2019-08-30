# Daniel Nichol and Jacob G Scott
# CA to simulate bacteria experiments.
# 15/4/2015
import numpy as np
import scipy as sp
from copy import copy
from random import randint
from random import random
from random import choice
from landscapes import oneStepNeighbours
from landscapes import unindex
from landscapes import convertGenotypeToInt

#Simulation constants.

p_divide = 1.0 #Probability to divide per time-step (time-step = 30mins, typical doubling time.)
N=4 #Genotype length
cc = 10**4 #carrying capacity
death_prob = 0.3 #probability of death per timestep.
mutation_prob = 10**-3
T_max = 10000
drug_conc = 0.3

# reproduction probability
def get_p_divide(MIC, drug_conc):
	prob = (1. - (drug_conc/MIC)) * p_divide
	return prob


def choose_random(pop):
	r = randint(1,sum(pop))
	s = 0
	gt = 0
	for i in range(len(pop)):
		s += pop[i]
		if r<=s:
			gt = i
			break

	return gt

landscape = [1.743, 1.662, 1.763, 1.785, 2.018, 2.05, 2.042, 0.218, 1.553, 0.256, 0.165, 0.221, 0.223, 0.239, 1.811, 0.288]
landscape = map(np.exp, landscape)
def update(pop, g_index, drug_conc):
	r1,r2,r3 = random(), random(), random()
	if pop[g_index] == 0:
		print "WTF IS GOING ON!?"
	delta_pop = [0 for i in range(16)]
	new_gi = g_index
	if r1 < death_prob:
		delta_pop[g_index] = -1
	elif sum(pop) < cc and r2<get_p_divide(landscape[g_index], drug_conc):
		if r3 < mutation_prob:
			osn = oneStepNeighbours(unindex(g_index,N))
			new_g = choice(osn)
			new_gi = convertGenotypeToInt(new_g)
	
		delta_pop[new_gi]+=1

	for k in range(len(pop)):
		pop[k] = pop[k] + delta_pop[k]

	return pop


def simpson_index(pop):
	tot = float(sum(pop))
	sq_props = map(lambda x : (x/tot)**2, pop)
	si = sum(sq_props)
	return 1-si

def shannon_entropy(pop):
	tot = float(sum(pop))
	shi = 0
	for p in pop:
		if p!=0:
			shi+=(p/tot)*np.log(p/tot)
	return -shi

import networkx as nx
import matplotlib.pyplot as plt

#Drawing machinary
G = nx.hypercube_graph(N)
pos=nx.spectral_layout(G) # positions for all nodes
labpos = {k : np.array([p[0],p[1]+0.1]) for k,p in pos.iteritems() }
base_node_size = 25.
node_scale = 1.5 #The scaling of node size by population at that node

#Recording information
population = [0 for i in range(2**N)]
population[0]=1000
life_history = [population]
div_history = [simpson_index(population)]
ent_history = [shannon_entropy(population)]

t = 0
while t < T_max:
	print "Timestep: ", t
	plt.clf()
	new_pop = copy(population)
	count_list = copy(population)
	while sum(count_list)!=0:
		g_index = choose_random(count_list)
		count_list[g_index]-=1
		new_pop = update(new_pop, g_index, drug_conc)
	life_history.append(population)
	div_history.append(simpson_index(population))
	ent_history.append(shannon_entropy(population))
	t+=1
	population = copy(new_pop)

	pop = life_history[-1]


	labels={}
	for i in range(len(pop)):
		labels[tuple(unindex(i,N))] = "".join(map(str,unindex(i,N)))
		node_color = 'b'
		if pop[i] > 0:
			node_color = 'r'

		plt.subplot(221)
		nx.draw_networkx_nodes(G, pos, nodelist=[tuple(unindex(i,N))], node_color=node_color, node_size=base_node_size+node_scale*pop[i], alpha=0.8)
		plt.axis('off')

		# plt.subplot(222)
		# nx.draw_networkx_nodes(G, pos, nodelist=[tuple(unindex(i,N))], node_color=node_color, node_size=base_node_size+node_scale*pop[i], alpha=0.8)
		# plt.axis('off')

		# plt.subplot(223)
		# nx.draw_networkx_nodes(G, pos, nodelist=[tuple(unindex(i,N))], node_color=node_color, node_size=base_node_size+node_scale*pop[i], alpha=0.8)
		# plt.axis('off')

		# plt.subplot(224)
		# nx.draw_networkx_nodes(G, pos, nodelist=[tuple(unindex(i,N))], node_color=node_color, node_size=base_node_size+node_scale*pop[i], alpha=0.8)
		# plt.axis('off')


	plt.subplot(221)
	nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.4)
	nx.draw_networkx_labels(G,labpos,labels,font_size=7)

	plt.subplot(222)
	plt.plot(div_history)
	plt.ylim(0.0,1.0)
	plt.xlim(0., T_max)

	plt.subplot(223)
	for i in range(2**N):
		plt.plot(np.array(life_history)[:,i])
	plt.ylim(0., cc)
	plt.xlim(0,T_max)


	plt.subplot(224)
	afs = []
	for pop in life_history:
		af = [0. for i in range(N)]
		for i in range(2**N):
			b = unindex(i,N)
			for j in range(len(b)):
				if b[j]==1:
					af[j]+=pop[i]
		af = map(lambda x : float(x)/sum(pop), af)
		afs.append(af)
	afs = np.array(afs).T
	for crv in afs:
		plt.plot(crv)
	plt.ylim(0,1.0)
	plt.xlim(0, T_max)


	# plt.subplot(222)
	# nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.4)
	# nx.draw_networkx_labels(G,labpos,labels,font_size=7)

	# plt.subplot(223)
	# nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.4)
	# nx.draw_networkx_labels(G,labpos,labels,font_size=7)

	# plt.subplot(224)
	# nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.4)
	# nx.draw_networkx_labels(G,labpos,labels,font_size=7)

	lab = str(t)
	lab = lab.zfill(4)
	plt.savefig("img"+lab+".png", dpi = 500)








#G.node[]
#nx.draw(G)
