from ctypes.wintypes import HACCEL
import random

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

class Conformation :

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(self.sequence)
        tab_temp1 = 2 * self.length * [0]
        tab_temp2 = self.length * [0]
        for i in range(len(tab_temp1)) :
            tab_temp1[i] = 2 * self.length * [0]
            tab_temp2[i // 2] = 2 * [0]
        self.representation = tab_temp1
        self.coordinate = tab_temp2

    # add covalent bond
    def printConformation(self) :
        for i in range(self.length * 2) :
            for j in range(self.length * 2) :
                if self.representation[i][j] != 0 :
                    print(self.representation[i][j], " ", end="")
                else :
                    print("   ", end="")
            print("")
        print("----------------------------------")

    def generateConformation(self) :
        i, j = self.length, self.length
        steps = [[-1,0], [1, 0], [0, -1], [0, 1]]
        for pos, residue in enumerate(self.sequence) :
            next_i, next_j = i, j
            while (self.representation[next_i][next_j] != 0) :
                step = random.choice(steps)
                next_i = i + step[0]
                next_j = j + step[1]
            i, j = next_i, next_j
            self.representation[i][j] = residue
            self.coordinate[pos] = [i, j]

    def translateToHP(self) :
        hydrophobic = ["A", "I", "L", "M", "F", "W", "Y", "V"]
        polar = ["T", "K", "R", "H"]
        sequence = list(self.sequence)
        for pos, residue in enumerate(sequence) :
            if residue in hydrophobic :
                sequence[pos] = "H"
            elif residue in polar :
                sequence[pos] = "P"
        self.sequence = sequence

    def __getFreePosition(self, k) :
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        i, j = self.coordinate[k][0], self.coordinate[k][1]
        for step in steps :
            if self.representation[i + step[0]][j + step[1]] == 0 :
                return [i + step[0], j + step[1]]
        return [0, 0]

    def endMoves (self, k) :
        if (k == 0 or k == self.length - 1) :
            i, j = self.coordinate[k][0], self.coordinate[k][1] 
            new_coordinate = self.__getFreePosition(k + 1) if k == 0 else self.__getFreePosition(k - 1) 
            new_i, new_j = new_coordinate[0], new_coordinate[1]
            self.representation[new_i][new_j] = self.representation[i][j]
            self.representation[i][j] = 0

    """

    def neighbourhoods () :
        return True


    def goodNeighbour(self, k, i, j) :
        previous = self.coordinate[k-1]
        next = self.coordinate[k+1]
        flag1 = False
        flag2 = False
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        for step in steps :
            if self.representation[i + step[0]][j + step[1]] == self.representation[previous[0]][previous[1]] :
                flag1 = True
            if self.representation[i + step[0]][j + step[1]] == self.representation[next[0]][next[1]] :
                flag2 = True
        return flag1 and flag2   
        
    def changeConformation(self, k) :
        i, j = self.coordinate[k][0], self.coordinate[k][1]
        next_i, next_j = i, j
        steps = [[-1,0], [1, 0], [0, -1], [0, 1]]
        flag = False
        for step in steps :
            next_i = i + step[0]
            next_j = j + step[1]
            if self.representation[next_i][next_j] == 0 and self.goodNeighbour(k, next_i, next_j) :
                flag = True
                break
        if flag : 
            self.representation[next_i][next_j] = self.representation[i][j]
            self.representation[i][j] = 0
        else :
            print("[Error] Unable to change residue positions : " + self.sequence[k])

    """ 


ma_conformation = Conformation("ARKLHGL")
#ma_conformation.translateToHP()
ma_conformation.generateConformation()
ma_conformation.printConformation()
ma_conformation.endMoves(6)
#ma_conformation.changeConformation(1)
ma_conformation.printConformation()
