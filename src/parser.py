"""Program that calculates the 2D/semi-3D folding of a protein using a REMC method.

Usage:
======
    python parser.py [-h] [-o OUTPUT] [-g GRAPHIC]
                 optimal n temp_min temp_max r p file [file ...]

    positionnal arguments :
        optimal: the optimal energy threshold
        n: the number of local steps in a Monte Carlo search
        temp_min: the minimum temperature values
        temp_max: the maximum temperature values
        r: the number of replicas to simulate
        file: a protein fasta file
    optional argument :
        -o OUTPUT, --output OUTPUT: directs the output to a name of your choice
        -g GRAPHIC, --graphic GRAPHIC: directs the output to a 2D or 3D graphic
"""

__authors__ = ("Moohmoo", "Université Paris Cité")
__copyright__ = "UPC"
__date__ = "2022-09-14"
__version__= "1.0.2"

import time
import argparse
import remc

def readFasta(name) :
    """Reads a fasta file and extracts the sequence.

    Parameters
    ----------
    name : string
        a parameter containing the name of the fasta file.

    Returns
    -------
    string
        return the sequence.
    """
    sequence = ""
    with open(name, "r") as filein :
        for line in filein :
            if (not line.startswith(">")) :
                sequence += line.rstrip()
    return sequence

parser = argparse.ArgumentParser(description='Calculates the 2D/semi-3D folding of a protein using a REMC method.')
parser.add_argument("optimal", type=int, help="the optimal energy threshold")
parser.add_argument("n", type=int, help="the number of local steps in a Monte Carlo search")
parser.add_argument("temp_min", type=int, help="the minimum temperature values")
parser.add_argument("temp_max", type=int, help="the maximum temperature values")
parser.add_argument("r", type=int, help="the number of replicas to simulate")
parser.add_argument("file", nargs='+', help="a protein fasta file")
parser.add_argument("-o", "--output", nargs=1, help="directs the output to a name of your choice")
parser.add_argument("-g", "--graphic", help="directs the output to a 2D or 3D graphic")
args = parser.parse_args()

for file_name in args.file :
    sequence = readFasta(file_name)
    start = time.time()
    random = True
    conformations = remc.REMCSimulation(sequence, args.optimal, args.n, args.temp_min, args.temp_max, args.r, random)
    best_conformation = remc.getBestConformation(conformations)
    print("Energy : ", best_conformation.getEnergy())
    end = time.time()
    print(f"The elapsed time in seconds : {end - start:0.2f}")
    args = vars(args)
    if (args["output"] is not None) :
        best_conformation.display2D(True, args["output"][0])
    if (args["graphic"] == "2d") :
        best_conformation.display2D(False, file_name+".png")
    elif (args["graphic"] == "3d") :
        best_conformation.display3D()

#python3 parser.py -g 2d -o config.png -10 500 220 220 2 ../test/P01013.fasta