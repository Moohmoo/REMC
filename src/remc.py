import random
import math
import copy
from tqdm import tqdm
import os

from conformation import Conformation

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

def MCsearch (n_mc, conformation) :
    print('process id:', os.getpid())
    for i in tqdm(range(n_mc), desc="MC research processing ") :
        temp_conformation = conformation.copy()
        k = random.uniform(0, temp_conformation.getLength() - 1)
        temp_conformation.changeConformation(k)
        temp_conformation.calculateEnergy()
        conformation.calculateEnergy()
        diff_energy = temp_conformation.getEnergy() - conformation.getEnergy()
        if diff_energy <= 0 :
            conformation = temp_conformation
        else :
            q = random.uniform(0, 1)
            if q > math.exp((-diff_energy)/conformation.getTemperature()) :
                conformation = temp_conformation
    return conformation

def REMCSimulation(sequence, optimal_energy, n_mc, temp_min, temp_max, n_replica, step, random) :
    temp_energy = 0
    conformations = [None for i in range(n_replica)]
    temperature = temp_min
    while temp_energy > optimal_energy :
        for i in range(n_replica) :
            conformation_temp = Conformation(sequence, temperature)
            temperature += ((temp_max - temp_min) // n_replica)
            if random :
                print("Generate a random conformation ...")
                conformation_temp.generateConformationRandom()
            else :
                print("Generate a linear conformation ...")
                conformation_temp.generateConformationLinear()
            conformations[i] = MCsearch(n_mc, conformation_temp)
            print("Energy conformation : ", conformations[i].getEnergy())
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
    
def main() :
    print('process id:', os.getpid())
    conformation = Conformation("ARKLHGLARKLHGLARKLHGLARKLHGL", 220)
    conformation.generateConformation()
    for i in tqdm(range(500), desc="MC research processing ") :
        temp_conformation = conformation.copy()
        k = random.uniform(0, temp_conformation.getLength() - 1)
        temp_conformation.changeConformation(k)
        temp_conformation.calculateEnergy()
        conformation.calculateEnergy()
        diff_energy = temp_conformation.getEnergy() - conformation.getEnergy()
        if diff_energy <= 0 :
            conformation = temp_conformation
        else :
            q = random.uniform(0, 1)
            if q > math.exp((-diff_energy)/conformation.getTemperature()) :
                conformation = temp_conformation
    print("Energy : ", conformation.getEnergy())
    conformation.generate3D()

if __name__ == "__main__" :
    main()