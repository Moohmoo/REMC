import random
import math

from residu import Residu

class Conformation :

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(self.sequence)
        self.__coordinate = self.length * [0, 0]
        len_representation = 2 * self.length + 1
        self.__representation = [[] for i in range(len_representation)]
        for i in range(len_representation) :
            self.__representation[i] = [Residu() for j in range(len_representation)]
        self.temperature = 15

    # add covalent bond
    def printConformation(self) :
        len_representation = 2 * self.length + 1
        for i in range(len_representation) :
            for j in range(len_representation) :
                print(self.__representation[i][j], " ", end="")
            print("")
        print("----------------------------------")

    # Prendre en compte les snail case
    def generateConformation(self) :
        i, j = self.length, self.length
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        previous, next = None, None
        for pos, residu in enumerate(self.sequence) :
            new_i, new_j = i, j
            while (self.__representation[new_i][new_j].getResidu() != None) :
                step = random.choice(steps)
                new_i = i + step[0]
                new_j = j + step[1]
            i, j = new_i, new_j
            self.__representation[i][j].setResidu(residu)
            self.__representation[i][j].setCoordinate([i, j])
            self.__coordinate[pos] = [i, j]
            
            if pos > 0 :
                next = self.__representation[i][j]
                previous.setNext(next)
                self.__representation[i][j].setPrevious(previous)
            previous = self.__representation[i][j]
            
    def translateToHP(self) :
        hydrophobicity = {"A": "H", "I": "H", "L": "H", "M": "H", "F": "H", 
        "W": "H", "Y": "H", "V": "H", "G": "H", "P": "H",
        "T": "P", "K": "P", "R": "P", "H": "P", "D": "P", "E": "P", "S": "P", 
        "N": "P", "Q": "P", "C": "P", "U": "P"}
        sequence = list(self.sequence)
        for k in range(len(sequence)) :
            i, j = self.__coordinate[k][0], self.__coordinate[k][1]
            self.__representation[i][j].setResidu(hydrophobicity[sequence[k]])

    def __getFreePosition(self, k) :
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        i, j = self.__coordinate[k][0], self.__coordinate[k][1]
        done = []
        while len(done) != 4 :
            step = random.choice(steps)
            if self.__representation[i + step[0]][j + step[1]].getResidu() == None :
                return [i + step[0], j + step[1]]
            if step not in done :
                done.append(step)
        return None

    def __move(self, new_coordinate, old_coordinate, k) :
        new_i, new_j = new_coordinate[0], new_coordinate[1]
        i, j = old_coordinate[0], old_coordinate[1]
        self.__representation[new_i][new_j].setResidu(self.__representation[i][j].getResidu())
        self.__representation[i][j].setResidu(None)
        self.__coordinate[k] = [new_i, new_j]
        print("Move done")

    def endMove(self, k) :
        if (k == 0 or k == self.length - 1) :
            #i, j = self.__coordinate[k][0], self.__coordinate[k][1] 
            new_coordinate = self.__getFreePosition(k + 1) if k == 0 else self.__getFreePosition(k - 1) 
            if (new_coordinate != None) :
                """
                new_i, new_j = new_coordinate[0], new_coordinate[1]
                self.__representation[new_i][new_j].setResidu(self.__representation[i][j].getResidu())
                self.__representation[i][j].setResidu(None)
                self.__coordinate[k] = [new_i, new_j]
                """
                self.__move(new_coordinate, self.__coordinate[k], k)
            else :
                print("[Error] EndMove impossible 1")
        else : 
            print("[Error] EndMove impossible 2")

    def cornerMove(self, k) :
        if (k > 0 and k < self.length - 1) :
            dist = math.dist(self.__coordinate[k-1], self.__coordinate[k+1])
            if (dist == math.sqrt(2)) :
                if (self.__coordinate[k + 1][0] == self.__coordinate[k][0] or
                    self.__coordinate[k - 1][1] == self.__coordinate[k][1]) :
                    new_coordinate = [self.__coordinate[k - 1][0], self.__coordinate[k + 1][1]]
                else :
                    new_coordinate = [self.__coordinate[k + 1][0], self.__coordinate[k - 1][1]]

                if (self.__representation[new_coordinate[0]][new_coordinate[1]].getResidu() == None) :
                    #i, j = self.__coordinate[k][0], self.__coordinate[k][1]
                    """
                    new_i, new_j = new_coordinate[0], new_coordinate[1]
                    self.__representation[new_i][new_j].setResidu(self.__representation[i][j].getResidu())
                    self.__representation[i][j].setResidu(None)
                    self.__coordinate[k] = [new_i, new_j]
                    """
                    self.__move(new_coordinate, self.__coordinate[k], k)
                else :    
                    print("[Error] CornerMoves impossible 1")
            else :
                print("[Error] CornerMoves impossible 2")
        else :
            print("[Error] CornerMoves impossible 3")

    def crankshaftMove(self, k) :
        if (k >= 1 and k <= self.length - 2) :
            dist1 = math.dist(self.__coordinate[k-1], self.__coordinate[k+1])
            dist2 = math.dist(self.__coordinate[k], self.__coordinate[k+2])
            dist3 = math.dist(self.__coordinate[k-1], self.__coordinate[k+2])
            if (dist1 == math.sqrt(2) and dist2 == math.sqrt(2) and dist3 == 1) :
                new_i1 = self.__coordinate[k-1][0] + (self.__coordinate[k-1][0] - self.__coordinate[k][0])
                new_j1 = self.__coordinate[k-1][1] + (self.__coordinate[k-1][1] - self.__coordinate[k][1])
                new_i2 = self.__coordinate[k+2][0] + (self.__coordinate[k+2][0] - self.__coordinate[k+1][0])
                new_j2 = self.__coordinate[k+2][1] + (self.__coordinate[k+2][1] - self.__coordinate[k+1][1])
                if (self.__representation[new_i1][new_j1].getResidu() == None and self.__representation[new_i2][new_j2].getResidu() == None) :
                    """
                    i1, j1 = self.__coordinate[k][0], self.__coordinate[k][1]
                    i2, j2 = self.__coordinate[k + 1][0], self.__coordinate[k + 1][1]
                    
                    self.__representation[new_i1][new_j1].setResidu(self.__representation[i1][j1].getResidu())
                    self.__representation[i1][j1].setResidu(None)
                    self.__coordinate[k] = [new_i1, new_j1]

                    self.__representation[new_i2][new_j2].setResidu(self.__representation[i2][j2].getResidu())
                    self.__representation[i2][j2].setResidu(None)
                    self.__coordinate[k+1] = [new_i2, new_j2]
                    """
                    new_coordinate1 = [new_i1, new_j1]
                    new_coordinate2 = [new_i2, new_j2]
                    self.__move(new_coordinate1, self.__coordinate[k], k)
                    self.__move(new_coordinate2, self.__coordinate[k + 1], k + 1)
                else :
                    print("[Error] CrankshaftMove impossible 1")
            else :
                print("[Error] CrankshaftMove impossible 2")
        else :
            print("[Error] CrankshaftMove impossible 3")

    def __isConnect(self, k) :
        if (k > 0 and k < self.length - 1) :
            dist = math.dist(self.__coordinate[k - 1], self.__coordinate[k + 1])
            return dist == 2

    def pullMove(self, k) :
        if (k >= 1 and k <= self.length - 2) :
            new_coordinate1 = self.__getFreePosition(k)
            new_coordinate2 = self.__getFreePosition(k + 1)
            if (new_coordinate1 != None and new_coordinate2 != None and new_coordinate1 != new_coordinate2) :
                old_coordinate1 = self.__coordinate[k]
                old_coordinate2 = self.__coordinate[k-1]
                self.__move(new_coordinate1, old_coordinate1)
                self.__move(new_coordinate2, old_coordinate2)
                i = 1
                while (not self.__isConnect(k - i)) :
                    print("is_connect ne fonctionne pas")
        else :
            print("[Error] pullMove impossible")

    # Loi Uniforme pour choisir un mouvement
    def changeConformation(self) :
        return True

    def calculateEnergy(self) :
        return 1

    def getSequence(self) :
        return self.sequence
    
    def getLength(self) :
        return self.length

    def getTemperature(self) :
        return self.temperature


if __name__ == "__main__" :
    ma_conformation = Conformation("ARKLHGL")
    ma_conformation.generateConformation()
    #ma_conformation.translateToHP()
    ma_conformation.printConformation()
    #ma_conformation.endMove(6)
    #ma_conformation.cornerMove(5)
    ma_conformation.crankshaftMove(2)
    ma_conformation.printConformation()
