import random

from residu import Residu

# Paramters : (φ, τ Min , τ max , χ, ρ) with φ the number of local steps in a Monte Carlo search, 
# τ min , and τ max , are the minimum and maximum temperature values respectively, χ is the
# number of replicas to simulate 

class Conformation :

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(self.sequence)
        self.coordinate = self.length * [0, 0]
        self.representation = [[] for i in range(2 * self.length)]
        for i in range(2 * self.length) :
            self.representation[i] = [Residu() for j in range(2 * self.length)]

    # add covalent bond
    def printConformation(self) :
        for i in range(self.length * 2) :
            for j in range(self.length * 2) :
                print(self.representation[i][j], " ", end="")
            print("")
        print("----------------------------------")

    def generateConformation(self) :
        i, j = self.length, self.length
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        previous, next = None, None
        for pos, residu in enumerate(self.sequence) :
            new_i, new_j = i, j
            while (self.representation[new_i][new_j].getResidu() != None) :
                step = random.choice(steps)
                new_i = i + step[0]
                new_j = j + step[1]
            i, j = new_i, new_j
            self.representation[i][j].setResidu(residu)
            self.representation[i][j].setCoordinate([i, j])
            self.coordinate[pos] = [i, j]
            
            if pos > 0 :
                next = self.representation[i][j]
                previous.setNext(next)
                self.representation[i][j].setPrevious(previous)
            previous = self.representation[i][j]
            
    def translateToHP(self) :
        hydrophobicity = {"A": "H", "I": "H", "L": "H", "M": "H", "F": "H", 
        "W": "H", "Y": "H", "V": "H", "G": "H", "P": "H",
        "T": "P", "K": "P", "R": "P", "H": "P", "D": "P", "E": "P", "S": "P", 
        "N": "P", "Q": "P", "C": "P", "U": "P"}
        sequence = list(self.sequence)
        for k in range(len(sequence)) :
            i, j = self.coordinate[k][0], self.coordinate[k][1]
            self.representation[i][j].setResidu(hydrophobicity[sequence[k]])

    def __getFreePosition(self, k) :
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        i, j = self.coordinate[k][0], self.coordinate[k][1]
        for step in steps :
            if self.representation[i + step[0]][j + step[1]].getResidu() == None :
                return [i + step[0], j + step[1]]
        return None

    def endMoves (self, k) :
        if (k == 0 or k == self.length - 1) :
            i, j = self.coordinate[k][0], self.coordinate[k][1] 
            new_coordinate = self.__getFreePosition(k + 1) if k == 0 else self.__getFreePosition(k - 1) 
            if (new_coordinate != None) :
                new_i, new_j = new_coordinate[0], new_coordinate[1]
                self.representation[new_i][new_j].setResidu(self.representation[i][j].getResidu())
                self.representation[i][j].setResidu(None)
                self.coordinate[k] = [new_i, new_j]
            else :
                print("[Error] Move impossible")

    def cornerMoves(self, k) :
        if (k != 0 and k != self.length - 1) :
            i, j = self.coordinate[k][0], self.coordinate[k][1]
            new_coordinate_left = self.__getFreePosition(k - 1)
            new_coordinate_right = self.__getFreePosition(k + 1)
            

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
        
    def changeConformation(self, k) :self.sequence = "".join(sequence)
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

if __name__ == "__main__" :
    ma_conformation = Conformation("ARKLHGL")

    ma_conformation.generateConformation()
    #ma_conformation.translateToHP()
    ma_conformation.printConformation()
    ma_conformation.endMoves(6)
    #ma_conformation.changeConformation(1)
    ma_conformation.printConformation()
