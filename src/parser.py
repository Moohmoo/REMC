import time
import argparse
from tqdm import tqdm

from residu import Residu
from conformation import Conformation, MovementError
import remc

def readFasta(name) :
    sequence = ""
    with open(name, "r") as filein :
        for line in filein :
            if (not line.startswith(">")) :
                sequence += line.rstrip()
    return sequence

def writeOutFile(name, conformation, modele) :
    with open(name, "w") as fileout :
        len_representation = 2 * conformation.length + 1
        fileout.write("----------------------------------\n")
        for i in range(len_representation) :
            fileout.write("|")
            for j in range(len_representation) :
                if modele :
                    if (conformation.getRepresentation()[i][j].getResidu() != None) :
                        fileout.write(conformation.getRepresentation()[i][j].getHydrophobicity())
                        fileout.write(" ")
                    else :
                        fileout.write("   ")
                else :
                    if (conformation.getRepresentation()[i][j].getResidu() != None) :
                        fileout.write(conformation.getRepresentation()[i][j].getResidu())
                        fileout.write(" ")
                    else :
                        fileout.write("   ")
                    fileout.write(" ")
            fileout.write("|\n")
        fileout.write("----------------------------------\n")

parser = argparse.ArgumentParser(description='Calculates the 2D folding of a protein using a REMC method.')
parser.add_argument("n", type=int, help="the number of local steps in a Monte Carlo search")
parser.add_argument("temp_min", type=int, help="the minimum temperature values")
parser.add_argument("temp_max", type=int, help="the maximum temperature values")
parser.add_argument("r", type=int, help="number of replicas to simulate")
parser.add_argument("p", type=int, help="the probability of performing a pull move")
parser.add_argument("file", nargs='+', help="a protein fasta file")
parser.add_argument("-o", "--output", help="directs the output to a name of your choice")
args = parser.parse_args()

for file_name in args.file :
    sequence = readFasta(file_name)
    start = time.time()
    conformations = remc.REMCSimulation(sequence, -2, args.n, args.temp_min, args.r, args.p)
    best_conformation = remc.getBestConformation(conformations)
    print("Energy : ", best_conformation.getEnergy())
    end = time.time()
    print(f"The elapsed time in seconds : {end - start:0.2f}")
    args = vars(parser.parse_args())
    args = vars(parser.parse_args())
    print(args)
    if (args["output"] is not None) :
        writeOutFile(args["output"], best_conformation, False)
