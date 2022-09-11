import random
import math
import copy
import time
from tqdm import tqdm
from multiprocessing import Process

from conformation import Conformation

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

def MCsearch (n_mc, conformation) :
    for i in tqdm(range(n_mc), desc="Monte Carlo research processing ") :
        temp_conformation = copy.deepcopy(conformation)
        k = random.uniform(0, temp_conformation.getLength() - 1)
        temp_conformation.changeConformation(k)
        temp_conformation.calculateEnergy()
        conformation.calculateEnergy()
        diff_energy = temp_conformation.getEnergy() - conformation.getEnergy()
        if diff_energy <= 0 :
            conformation = copy.deepcopy(temp_conformation)
        else :
            q = random.uniform(0, 1)
            if q > math.exp((-diff_energy)/conformation.getTemperature()) :
                conformation = copy.deepcopy(temp_conformation) 
    return conformation

def REMCSimulation(sequence, optimal_energy, n_mc, temperature, n_replica, step) :
    temp_energy = 0
    conformations = [None for i in range(n_replica)]
    while temp_energy > optimal_energy :
        for i in range(n_replica) :
            conformation_temp = Conformation(sequence, temperature)
            conformation_temp.generateConformation()
            p = Process(target=MCsearch, args=(n_mc, conformation_temp,))
            p.start()
            p.join()
            conformations[i] = MCsearch(n_mc, conformation_temp)
            if conformations[i].getEnergy() < temp_energy :
                temp_energy = conformations[i].getEnergy()
        
        for i in range(n_replica - 1) :
            j = i + 1
            betaj = 1/conformations[j].getTemperature()
            betai = 1/conformations[i].getTemperature()
            energyi = conformations[i].getEnergy()
            energyj = conformations[j].getEnergy()
            delta = (betaj - betai)*(energyi - energyj)
            if delta <= 0 :
                temperature_temp = conformations[i].getTemperature()
                conformations[i].setTemperature(conformations[j])
                conformations[j].setTemperature(temperature_temp)
            else :
                q = random.uniform(0, 1)
                if q < math.exp(-delta) : 
                    temperature_temp = conformations[i].getTemperature()
                    conformations[i].setTemperature(conformations[j])
                    conformations[j].setTemperature(temperature_temp)
    print("[REMC Simulation done]")
    return conformations
    
def getBestConformation(conformations) :
    best_conformation = conformations[0]
    for i in range(1, len(conformations)) :
        if (best_conformation.getEnergy() > conformations[i].getEnergy()) :
            best_conformation = conformations[i]
    return best_conformation

if __name__ == "__main__" :
    start = time.time()
    conformations = REMCSimulation("ARKLHGLARKLHGLARKLHGLARKLHGL", -2, 500, 220, 2, 0.5)
    best_conformation = getBestConformation(conformations)
    best_conformation.printConformation(False)
    print("Energy : ", best_conformation.getEnergy())
    end = time.time()
    print(f"The elapsed time in seconds : {end - start:0.2f}")