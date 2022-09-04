import random

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

class Conformation :

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(self.sequence)
        tab_temp1 = 2*self.length*[0]
        tab_temp2 = self.length*[0]
        for i in range(len(tab_temp1)) :
            tab_temp1[i] = 2*self.length*[0]
            tab_temp2[i//2] = 2*[0]
        self.representation = tab_temp1
        self.coordinate = tab_temp2

    def generateConformation(self) :
        i, j = self.length, self.length
        for pos, residue in enumerate(self.sequence) :
            while (self.representation[i][j] != 0) :
                i = i + random.randint(-1, 1)
                j = j + random.randint(-1, 1)
            self.representation[i][j] = residue
            self.coordinate[pos] = [i, j]

    def printConformation(self) :
        for i in range(self.length * 2) :
            for j in range(self.length * 2) :
                print(self.representation[i][j], " ", end="")
            print("")
        print("----------------------------------")

    def isAdjacent(self, i, j) :
        if (self.representation[i+1][j] or self.representation[i-1][j] or self.representation[i][j+1] or self.representation[i][j-1]
        or self.representation[i+1][j+1] or self.representation[i+1][j-1] or self.representation[i-1][j+1] or self.representation[i-1][j-1]) :
            return True
        else :
            return False
    
    def isAllAdjacent(self) :
        for k in range(len(self.sequence)) :
            i, j = self.coordinate[k][0], self.coordinate[k][1]
            if (not self.isAdjacent(i, j)) :
                return False
        return True

    def changeConformation(self, k) :
        temp = Conformation(self.sequence)
        temp.representation = self.representation
        temp.coordinate = self.coordinate

        previousI, previousJ = self.coordinate[k-1][0], self.coordinate[k-1][1] 
        i, j = previousI, previousJ
        while (True) :
            print("test")
            while (self.representation[i][j] != 0 or (not self.isAdjacent(i, j))) :
                i = i + random.randint(-1, 1)
                j = j + random.randint(-1, 1)
        
            temp.representation[i][j] = temp.representation[previousI][previousJ]
            temp.representation[previousI][previousJ] = 0

            if temp.isAllAdjacent() :
                break

        self.representation[i][j] = self.representation[previousI][previousJ]
        self.representation[previousI][previousJ] = 0

        


ma_conformation = Conformation("ARKLHGL")
ma_conformation.generateConformation()
ma_conformation.printConformation()
ma_conformation.changeConformation(2)
ma_conformation.printConformation()
