"""Python file grouping the functions necessary for the REMC (replica exchange Monte-Carlo) algorithm
"""

from conformation import Conformation

import random
import math
import copy
from tqdm import tqdm
import os

def MCsearch (n_mc, conformation) :
    """Performs the Monte-Carlo search explained in the article (by Thachuk, C. and al.).

    Parameters
    ----------
    n_mc : int
        the number of local steps in a Monte Carlo search.

    Returns
    -------
    object
        return the object corresponding to a conformation after a simulation.
    """
    print('process id:', os.getpid())
    for i in tqdm(range(n_mc), desc="MC research processing ") :
        temp_conformation = copy.deepcopy(conformation)
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

def REMCSimulation(sequence, optimal_energy, n_mc, temp_min, temp_max, n_replica, random) :
    """Performs the replica exchange Monte-Carlo explained in the article (by Thachuk, C. and al.).

    Parameters
    ----------
    sequence : str
        the string containing the sequence.
    optimal energy : int
        the optimal energy threshold.
    n_mc : int
        the number of local steps in a Monte Carlo search.
    temp_min : int
        the minimum temperature values.
    temp_max : int
        the maximum temperature values.
    n_replica : int
        the number of replicas to simulate.
    random : boolean
        the type of conformation to generate.

    Returns
    -------
    list
        return a list containing the replicas after the simulations.
    """
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
            print("Energy of replica ", i, " : ", conformations[i].getEnergy())
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
    print("REMC Simulation done !")
    return conformations
    
def getBestConformation(conformations) :
    """Determines the best conformation among a list of conformations.

    Parameters
    ----------
    conformations : list
        a list containing conformation objects.

    Returns
    -------
    object
        return the object corresponding to the best conformation.
    """
    best_conformation = conformations[0]
    for i in range(1, len(conformations)) :
        if (best_conformation.getEnergy() > conformations[i].getEnergy()) :
            best_conformation = conformations[i]
    return best_conformation
