import random
import math
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


from residu import Residu

class Conformation :

    def __init__(self, sequence, temperature):
        self.sequence = sequence
        self.length = len(self.sequence)
        self.__coordinate = [[0, 0] for i in range(self.length)]
        self._createRepresentation()
        self.temperature = temperature
        self.energy = 0

    def _createRepresentation(self) :
        len_representation = 2 * self.length + 1
        self.__representation = [[] for i in range(len_representation)]
        for i in range(len_representation) :
            self.__representation[i] = [Residu() for j in range(len_representation)]

    # add covalent bond
    def printConformation(self, modele) :
        len_representation = 2 * self.length + 1
        print("----------------------------------")
        for i in range(len_representation) :
            print("|", end="")
            for j in range(len_representation) :
                if modele :
                    if (self.__representation[i][j].getResidu() != None) :
                        print(self.__representation[i][j].getHydrophobicity(), " ", end="")
                    else :
                        print("   ", end="")
                else :
                    print(self.__representation[i][j], " ", end="")
            print("|")
        print("----------------------------------")

    def generate3D(self) :
        fig = plt.figure(figsize=(5,5))
        ax = plt.axes(projection='3d')
        plt.grid(False)
        plt.axis(False)
        coordinate = self.__coordinate
        x, y = [], []
        previousX = coordinate[0][0]
        previousY = coordinate[0][1]
        for i in range(1, len(coordinate)) :
            x.append(coordinate[i][0])
            y.append(coordinate[i][1])
            ax.plot([0, 0] , [previousX, coordinate[i][0]], [previousY, coordinate[i][1]], color="black")
            ax.text(0, coordinate[i][0] + 0.15, coordinate[i][1] + 0.25, "H")
            previousX = coordinate[i][0]
            previousY = coordinate[i][1]
        ax.scatter(0, x, y)
        plt.show()

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

    # Prendre en compte les snail case
    def generateConformation(self) :
        i, j = self.length, self.length
        self.__representation[i][j].setResidu(self.sequence[0])
        self.__coordinate[0] = [i, j]

        for k, residu in enumerate(self.sequence[1:], 1) :
            new_coordinate = self.__getFreePosition(k - 1)
            if (new_coordinate == None) :
                print("Snail case detected")
                self._createRepresentation()
                self.generateConformation()
                return
            else :
                i, j = new_coordinate[0], new_coordinate[1]
                self.__representation[i][j].setResidu(residu)
                self.__representation[i][j].translateHP()
                self.__coordinate[k] = [i, j]
        """
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        previous, next = None, None

        for k, residu in enumerate(self.sequence) :
            new_i, new_j = i, j
            while (self.__representation[new_i][new_j].getResidu() != None) :
                step = random.choice(steps)
                new_i = i + step[0]
                new_j = j + step[1]
            i, j = new_i, new_j
            self.__representation[i][j].setResidu(residu)
            self.__representation[i][j].setCoordinate([i, j])
            self.__coordinate[k] = [i, j]
        """
        """
            if k > 0 :
                next = self.__representation[i][j]
                previous.setNext(next)
                self.__representation[i][j].setPrevious(previous)
            previous = self.__representation[i][j]
        """

    def checkCoordinate(self) :
        for k in range(len(self.__coordinate)) :
            i, j = self.__coordinate[k][0], self.__coordinate[k][1]
            if i >= (2 * self.length) or j >= (2 * self.length) or i == 0 or j == 0 :
                self._centerPosition(i, j)
                return

    def __move(self, new_coordinate, old_coordinate, k) :
        new_i, new_j = new_coordinate[0], new_coordinate[1]
        i, j = old_coordinate[0], old_coordinate[1]
        self.__representation[new_i][new_j].setResidu(self.__representation[i][j].getResidu())
        self.__representation[new_i][new_j].translateHP()
        self.__representation[i][j].setResidu(None)
        self.__representation[i][j].setHydrophobicity(None)
        self.__coordinate[k] = [new_i, new_j]

    # A MODIFIER
    def _centerPosition(self, i, j) :
        if (i >= (2 * self.length)) :
            for k in range(len(self.__coordinate)) :
                i, j = self.__coordinate[k][0], self.__coordinate[k][1]
                new_i, new_j = i, j

                new_i =  i - self.length

                self.__move([new_i, new_j], self.__coordinate[k], k)
        elif (j >= (2 * self.length)) :
            for k in range(len(self.__coordinate)) :
                i, j = self.__coordinate[k][0], self.__coordinate[k][1]
                new_i, new_j = i, j

                new_j =  j - self.length

                self.__move([new_i, new_j], self.__coordinate[k], k)
        elif (i == 0) :
            for k in range(len(self.__coordinate)) :
                i, j = self.__coordinate[k][0], self.__coordinate[k][1]
                new_i, new_j = i, j

                new_i =  i + self.length

                self.__move([new_i, new_j], self.__coordinate[k], k)
        else :
            for k in range(len(self.__coordinate)) :
                i, j = self.__coordinate[k][0], self.__coordinate[k][1]
                new_i, new_j = i, j

                new_j =  j + self.length
                self.__move([new_i, new_j], self.__coordinate[k], k)

    def endMove(self, k) :
        if (k == 0 or k == self.length - 1) :
            new_coordinate = self.__getFreePosition(k + 1) if k == 0 else self.__getFreePosition(k - 1) 
            if (new_coordinate != None) :
                self.__move(new_coordinate, self.__coordinate[k], k)
            else :
                raise MovementError("EndMove impossible")
        else :
            raise MovementError("EndMove impossible")

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
                    self.__move(new_coordinate, self.__coordinate[k], k)
                else :
                    raise MovementError("CornerMoves impossible")
            else :
                raise MovementError("CornerMoves impossible")
        else :
            raise MovementError("CornerMoves impossible")

    def crankshaftMove(self, k) :
        if (k >= 1 and k < self.length - 2) :
            dist = math.dist(self.__coordinate[k-1], self.__coordinate[k+2])
            if (dist == 1) :
                new_i1 = self.__coordinate[k-1][0] + (self.__coordinate[k-1][0] - self.__coordinate[k][0])
                new_j1 = self.__coordinate[k-1][1] + (self.__coordinate[k-1][1] - self.__coordinate[k][1])
                new_i2 = self.__coordinate[k+2][0] + (self.__coordinate[k+2][0] - self.__coordinate[k+1][0])
                new_j2 = self.__coordinate[k+2][1] + (self.__coordinate[k+2][1] - self.__coordinate[k+1][1])
                if (self.__representation[new_i1][new_j1].getResidu() == None and self.__representation[new_i2][new_j2].getResidu() == None) :
                    new_coordinate1 = [new_i1, new_j1]
                    new_coordinate2 = [new_i2, new_j2]
                    self.__move(new_coordinate1, self.__coordinate[k], k)
                    self.__move(new_coordinate2, self.__coordinate[k + 1], k + 1)
                else :
                    raise MovementError("CrankshaftMove impossible")
            else :
                raise MovementError("CrankshaftMove impossible")
        else :
            raise MovementError("CrankshaftMove impossible")

    def __isConnect(self, k) :
        if (k > 0 and k < self.length - 1) :
            dist1 = math.dist(self.__coordinate[k], self.__coordinate[k - 1])
            dist2 = math.dist(self.__coordinate[k], self.__coordinate[k + 1])
            return dist1 == 1 and dist2 == 1
        return None

    def pullMove(self, k) :
        sens = random.choice([0, 1])
        positions = [[k + 1, k - 1, 2, 1], [k - 1, k + 1, -2, -1]]
        if (k >= 1 and k <= self.length - 2) :
            new_coordinate1 = self.__getFreePosition(k)
            #new_coordinate2 = self.__getFreePosition(k + 1)
            new_coordinate2 = self.__getFreePosition(positions[sens][0])
            if (new_coordinate1 != None and new_coordinate2 != None and new_coordinate1 != new_coordinate2 and 
                math.dist(new_coordinate1, new_coordinate2) == 1 and 
                math.dist(new_coordinate1, self.__coordinate[positions[sens][1]]) == math.sqrt(2) and
                math.dist(new_coordinate2, self.__coordinate[k]) == math.sqrt(2)) :
                old_coordinate1 = self.__coordinate[positions[sens][1]]
                old_coordinate2 = self.__coordinate[k]
                self.__move(new_coordinate2, old_coordinate2, k)
                self.__move(new_coordinate1, old_coordinate1, positions[sens][1])
                i = positions[sens][2]
                while (k - i >= 0 and not self.__isConnect(k - i)) :
                    self.__move(old_coordinate2, self.__coordinate[k - i], k - i)
                    old_coordinate2 = old_coordinate1
                    i += positions[sens][3]
            else :
                raise MovementError("PullMove impossible")
        else :
            raise MovementError("PullMove impossible")

    def changeConformation(self, k) :
        try :
            if k == 0 or k == self.length - 1 :
                self.endMove(k)
            else :
                case = int(random.uniform(0, 2))
                if case == 0 :
                    self.cornerMove(k)
                elif case == 1 :
                    self.crankshaftMove(k)
                else :
                    self.pullMove(k)
        except :
            #print("The conformation has not been changed")
            """
            k = random.choice(list(range(0, self.length)))
            self.changeConformation(k)
            """
            k = random.choice(list(range(0, self.length)))
            self.changeConformation(k)
            
    def calculateEnergy(self) :
        energy = 0
        for k in range(len(self.__coordinate)) :
            for l in range(k, len(self.__coordinate)) :
                ik, jk = self.__coordinate[k][0], self.__coordinate[k][1]
                il, jl = self.__coordinate[l][0], self.__coordinate[l][1]
                dist = math.dist(self.__coordinate[k], self.__coordinate[l])
                hydrophobick = self.__representation[ik][jk].getHydrophobicity()
                hydrophobicl = self.__representation[il][jl].getHydrophobicity()
                if (dist == 1 and abs(k - l) != 1 and hydrophobick == hydrophobicl and hydrophobick == "H") :
                    energy -= 1
        self.energy = energy

    def getResidu(self, k) :
        i, j = self.__coordinate[k][0], self.__coordinate[k][1]
        return self.__representation[i][j].getResidu()

    def getHydrophicity(self, k) :
        i, j = self.__coordinate[k][0], self.__coordinate[k][1]

    def getSequence(self) :
        return self.sequence
    
    def getLength(self) :
        return self.length

    def getTemperature(self) :
        return self.temperature
    
    def getEnergy(self) :
        return self.energy

    def setTemperature(self, temperature) :
        self.temperature = temperature
        
    # A supprimer
    def getCoordinate(self) :
        return self.__coordinate

    def getRepresentation(self) :
        return self.__representation

class MovementError(Exception):
    """Raised when a movement is impossible"""
    pass
"""
ma_conformation = Conformation("ARKLHGLARKLHGLARKLHGL", 220)
ma_conformation.generateConformation()
ma_conformation.printConformation(False)
for i in tqdm(range(500), "TEST MC ") :
#for i in range(500) :
    k = random.uniform(0, ma_conformation.getLength() - 1)
    #temp_conformation = ma_conformation
    #temp_conformation = copy.deepcopy(ma_conformation)
    temp_conformation = ujson.loads(ujson.dumps(ma_conformation))
    temp_conformation.changeConformation(k)
    temp_conformation.calculateEnergy()
    ma_conformation.calculateEnergy()
    #print(ma_conformation.getEnergy(), " ", temp_conformation.getEnergy())
ma_conformation.printConformation(False)
"""
"""
if __name__ == "__main__" :
    modele = False
    ma_conformation = Conformation("ARKLHGLARKLHGLARKLHGLARKLHGLAR", 220)
    ma_conformation.generateConformation()
    ma_conformation.translateHP()
    ma_conformation.printConformation(modele)
    print("Len sequence : ", ma_conformation.length, " | Taille representation : ", len(ma_conformation.getRepresentation()[0]))
    print("Done !")
    ma_conformation.translateHP()
    ma_conformation.printConformation(modele)
    for i in range (10000) :
        k = random.choice(list(range(0, ma_conformation.length)))
        ma_conformation.changeConformation(k)
        ma_conformation.checkCoordinate()
    ma_conformation.printConformation(False)
    ma_conformation.printConformation(True)
    ma_conformation.calculateEnergy()
    ma_conformation.getEnergy()
    print(ma_conformation.length, " ", len(ma_conformation.getRepresentation[0]))
    print("Done !")
    #ma_conformation.calculateEnergy()
    #ma_conformation.endMove(6)
    #ma_conformation.cornerMove(5)
    #ma_conformation.crankshaftMove(2)
    #ma_conformation.pullMove(3)
    #ma_conformation.printConformation()
"""
if __name__ == "__main__" :
    modele = False
    ma_conformation = Conformation("ARKLHGLARKLHGLARKLHGLARKLHGLAR", 220)
    ma_conformation.generateConformation()
    ma_conformation.generate3D()
