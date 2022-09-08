import random
import math
import copy

from residu import Residu
from conformation import Conformation

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

def MCsearch (n_mc, conformation) :
    for i in range(n_mc) :
        temp_conformation = copy.copy(conformation)
        k = random.uniform(0, temp_conformation.getLength() - 1)
        temp_conformation.changeConformation()
        diff_energy = temp_conformation.calculateEnergy() - conformation.calculateEnergy()
        if diff_energy <= 0 :
            conformation = copy.copy(temp_conformation)
        else :
            q = random.uniform(0, 1)
            if q > math.exp((-diff_energy)/conformation.getTemprature()) :
                conformation = copy.copy(temp_conformation) 
    return conformation

def REMCSimulation(n_replica, n_mc) :
    optimal_energy = 0
    temp_energy = 0
    conformations = [[] for i in range(n_replica)]
    while temp_energy > optimal_energy :
        for i in range(n_replica) :
            conformations[i].append(MCsearch(n_mc))
    return True