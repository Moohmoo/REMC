"""Python file grouping the Conformation and MovementError classes as well as the associated functions.
"""

from residu import Residu

import random
import math
import matplotlib.pyplot as plt

class Conformation :
    """A class representing a 2D conformation of a given protein.

    Attributes
    ----------
        sequence : string
            protein sequence
        length : int
            the length of the sequence
        coordinate : list
            list containing the coordinates of each residu
        residus : list
            list containing residu objects
        temperature : int
            a temperature of the protein
        energy : int
            protein energy
    """

    def __init__(self, sequence, temperature):
        """The constructor of the Conformation object

        Parameters
        ----------
        self : object
            a parameter in instance method.
        sequence : string
            a parameter corresponding to a sequence of a protein
        temperature : string
            a parameter corresponding to a temperature of a protein

        Returns
        -------
        None
            return nothing.
        """
        self.sequence = sequence
        self.length = len(self.sequence)
        self.__coordinate = [[0, 0] for i in range(self.length)]
        self.residus = [Residu() for i in range(self.length)]
        self.temperature = temperature
        self.energy = 0

    def display2D(self, save, name) :
        """Displays a 2d representation of conformation.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        save : boolean
            a parameter allowing to know if the figure must be saved.
        name : str
            parameter storing the name of the file to save.

        Returns
        -------
        None
            return nothing.
        """
        fig = plt.figure(figsize=(5,5))
        ax = plt.axes()
        plt.grid(False)
        plt.axis(False)
        coordinate = self.__coordinate
        xh, yh = [], []
        xp, yp = [], []
        if self.residus[0].getHydrophobic() == "H" :
            xh.append(coordinate[0][0])
            yh.append(coordinate[0][1])
        else :
            xp.append(coordinate[0][0])
            yp.append(coordinate[0][1])
        ax.text(coordinate[0][0] + 0.1, coordinate[0][1] + 0.15, self.residus[0].getResidu())
        previousX = coordinate[0][0]
        previousY = coordinate[0][1]
        for i in range(1, len(coordinate)) :
            plt.plot([previousX, coordinate[i][0]], [previousY, coordinate[i][1]], color="black")
            ax.text(coordinate[i][0] + 0.1, coordinate[i][1] + 0.15, self.residus[i].getResidu())
            previousX = coordinate[i][0]
            previousY = coordinate[i][1]
            if (self.residus[i].getHydrophobic() == "H") :
                xh.append(coordinate[i][0])
                yh.append(coordinate[i][1])
            else :
                xp.append(coordinate[i][0])
                yp.append(coordinate[i][1])
        plt.scatter(xh, yh, color="red", s=50, label="Hydrophobic", zorder=2)
        plt.scatter(xp, yp, color="blue", s=50, label="Polar", zorder=2)
        plt.legend()
        if (save) :
            directory = "../res/"+name
            fig.savefig(directory)
        else :
            plt.show()
        plt.close()

    def display3D(self) :
        """Displays a 3d representation of conformation.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        fig = plt.figure(figsize=(5,5))
        ax = plt.axes(projection='3d')
        plt.grid(False)
        plt.axis(False)
        coordinate = self.__coordinate
        xh, yh = [], []
        xp, yp = [], []
        if self.residus[0].getHydrophobic() == "H" :
            xh.append(coordinate[0][0])
            yh.append(coordinate[0][1])
        else :
            xp.append(coordinate[0][0])
            yp.append(coordinate[0][1])
        ax.text(0, coordinate[0][0] + 0.15, coordinate[0][1] + 0.25, self.residus[0].getResidu())
        previousX = coordinate[0][0]
        previousY = coordinate[0][1]
        for i in range(1, len(coordinate)) :
            ax.plot([0, 0] , [previousX, coordinate[i][0]], [previousY, coordinate[i][1]], color="black")
            ax.text(0, coordinate[i][0] + 0.15, coordinate[i][1] + 0.25, self.residus[i].getResidu())
            previousX = coordinate[i][0]
            previousY = coordinate[i][1]
            if (self.residus[i].getHydrophobic() == "H") :
                xh.append(coordinate[i][0])
                yh.append(coordinate[i][1])
            else :
                xp.append(coordinate[i][0])
                yp.append(coordinate[i][1])
        ax.scatter(0, xh, yh, color="red", s=50, label="Hydrophobic", zorder=2)
        ax.scatter(0, xp, yp, color="blue", s=50, label="Polar", zorder=2)
        plt.legend()
        plt.show()
        plt.close()

    def __getFreePosition(self, k) :
        """Determines the available adjacent positions of a residue k.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu

        Returns
        -------
        list/None
            returns a list containing the coordinates found or None.
        """
        steps = [[-1, 0], [1, 0], [0, -1], [0, 1]]
        i, j = self.__coordinate[k][0], self.__coordinate[k][1]
        done = []
        while len(done) != 4 :
            step = random.choice(steps)
            new_coordinate = [i + step[0], j + step[1]]
            if new_coordinate not in self.__coordinate :
                return new_coordinate
            if step not in done :
                done.append(step)
        return None

    def generateConformationLinear(self) :
        """Generates the positions of the residues randomly.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        for k in range(self.length) :
            i, j = self.length, self.length
            self.residus[k].setResidu(self.sequence[k])
            self.residus[k].translateHP()
            self.__coordinate[k] = [i, j + k]

    def generateConformationRandom(self) :
        """Generates the positions of the residues randomly.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        i, j = self.length, self.length
        self.residus[0].setResidu(self.sequence[0])
        self.residus[0].translateHP()
        self.__coordinate[0] = [i, j]

        for k, residu in enumerate(self.sequence[1:], 1) :
            new_coordinate = self.__getFreePosition(k - 1)
            if (new_coordinate == None) :
                self.__coordinate = [[0, 0] for i in range(self.length)]
                self.generateConformationRandom()
                return
            else :
                i, j = new_coordinate[0], new_coordinate[1]
                self.residus[k].setResidu(residu)
                self.residus[k].translateHP()
                self.__coordinate[k] = [i, j]
       
    def __move(self, new_coordinate, k) :
        """Changes the position of the residue k.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        new_coordinate : list
            a parameter that contains the new coordinates.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.
        """
        self.__coordinate[k] = [new_coordinate[0], new_coordinate[1]]

    def endMove(self, k) :
        """Makes a movement of a residue from the beginning or the end.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.

        Raises
        ------
        MovementError
            raises an exception with a message.
        """
        if (k == 0 or k == self.length - 1) :
            new_coordinate = self.__getFreePosition(k + 1) if k == 0 else self.__getFreePosition(k - 1) 
            if (new_coordinate != None) :
                self.__move(new_coordinate, k)
            else :
                raise MovementError("EndMove impossible")
        else :
            raise MovementError("EndMove impossible")

    def cornerMove(self, k) :
        """Makes a corner move of a random residue excluding start and end.
        
        The movement is explained on the article (by Thachuk, C. and al.).

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.

        Raises
        ------
        MovementError
            raises an exception with a message.
        """
        if (k > 0 and k < self.length - 1) :
            dist = math.dist(self.__coordinate[k-1], self.__coordinate[k+1])
            if (dist == math.sqrt(2)) :
                if (self.__coordinate[k + 1][0] == self.__coordinate[k][0] or
                    self.__coordinate[k - 1][1] == self.__coordinate[k][1]) :
                    new_coordinate = [self.__coordinate[k - 1][0], self.__coordinate[k + 1][1]]
                else :
                    new_coordinate = [self.__coordinate[k + 1][0], self.__coordinate[k - 1][1]]

                if new_coordinate not in self.__coordinate :
                    self.__move(new_coordinate, k)
                else :
                    raise MovementError("CornerMoves impossible")
            else :
                raise MovementError("CornerMoves impossible")
        else :
            raise MovementError("CornerMoves impossible")

    def crankshaftMove(self, k) :
        """Makes a crankshaft move of a random residue excluding start and end.
        
        The movement is explained on the article (by Thachuk, C. and al.).

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.

        Raises
        ------
        MovementError
            raises an exception with a message.
        """
        if (k >= 1 and k < self.length - 2) :
            dist = math.dist(self.__coordinate[k-1], self.__coordinate[k+2])
            if (dist == 1) :
                new_i1 = self.__coordinate[k-1][0] + (self.__coordinate[k-1][0] - self.__coordinate[k][0])
                new_j1 = self.__coordinate[k-1][1] + (self.__coordinate[k-1][1] - self.__coordinate[k][1])
                new_i2 = self.__coordinate[k+2][0] + (self.__coordinate[k+2][0] - self.__coordinate[k+1][0])
                new_j2 = self.__coordinate[k+2][1] + (self.__coordinate[k+2][1] - self.__coordinate[k+1][1])
                if ([new_i1, new_j1] not in self.__coordinate and [new_i2, new_j2] not in self.__coordinate) :
                    new_coordinate1 = [new_i1, new_j1]
                    new_coordinate2 = [new_i2, new_j2]
                    self.__move(new_coordinate1, k)
                    self.__move(new_coordinate2, k + 1)
                else :
                    raise MovementError("CrankshaftMove impossible")
            else :
                raise MovementError("CrankshaftMove impossible")
        else :
            raise MovementError("CrankshaftMove impossible")

    def __isConnect(self, k) :
        """Checks if two residus are side by side.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        boolean/None
            return True if distance between 2 residu equal 1 else nothing.
        """
        if (k > 0 and k < self.length - 1) :
            dist1 = math.dist(self.__coordinate[k], self.__coordinate[k - 1])
            dist2 = math.dist(self.__coordinate[k], self.__coordinate[k + 1])
            return dist1 == 1 and dist2 == 1
        return None

    def pullMove(self, k) :
        """Makes a pull move of a random residue excluding start and end.
        
        The movement is explained on the article (by Thachuk, C. and al.).

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.

        Raises
        ------
        MovementError
            raises an exception with a message.
        """
        sens = random.choice([0, 1])
        positions = [[k + 1, k - 1, 2, 1], [k - 1, k + 1, -2, -1]]
        if (k >= 1 and k <= self.length - 2) :
            new_coordinate1 = self.__getFreePosition(k)
            new_coordinate2 = self.__getFreePosition(positions[sens][0])
            if (new_coordinate1 != None and new_coordinate2 != None and new_coordinate1 != new_coordinate2 and 
                math.dist(new_coordinate1, new_coordinate2) == 1 and 
                math.dist(new_coordinate1, self.__coordinate[positions[sens][1]]) == math.sqrt(2) and
                math.dist(new_coordinate2, self.__coordinate[k]) == math.sqrt(2)) :
                old_coordinate1 = self.__coordinate[positions[sens][1]]
                old_coordinate2 = self.__coordinate[k]
                self.__move(new_coordinate2, k)
                self.__move(new_coordinate1, positions[sens][1])
                i = positions[sens][2]
                while (k - i >= 0 and not self.__isConnect(k - i)) :
                    self.__move(old_coordinate2, k - i)
                    old_coordinate2 = old_coordinate1
                    i += positions[sens][3]
            else :
                raise MovementError("PullMove impossible")
        else :
            raise MovementError("PullMove impossible")

    def changeConformation(self, k) :
        """Applies the different movements on the conformation.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        None
            return nothing.
        """
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
            k = random.choice(list(range(0, self.length)))
            self.changeConformation(k)
            
    def calculateEnergy(self) :
        """Calculate the energy of a conformation according to the article.

        The calculation is explained on the article (by Thachuk, C. and al.).

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        energy = 0
        for k in range(len(self.__coordinate)) :
            for l in range(k, len(self.__coordinate)) :
                dist = math.dist(self.__coordinate[k], self.__coordinate[l])
                hydrophobick = self.residus[k].getHydrophobic()
                hydrophobicl = self.residus[l].getHydrophobic()
                if (dist == 1 and abs(k - l) != 1 and hydrophobick == hydrophobicl and hydrophobick == "H") :
                    energy -= 1
        self.energy = energy

    def getResidu(self, k) :
        """Determines the residu according to a given position.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        str
            return a string corresponding to the residu.
        """
        return self.residus[k].getResidu()

    def getHydrophicity(self, k) :
        """Determines the polarity of a residu according to a given position.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        k : int
            index of a residu.

        Returns
        -------
        str
            return a string corresponding to the polarity.
        """
        return self.residus[k].getHydrophobic()

    def getSequence(self) :
        """Determines the sequence.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        str
            return a string corresponding to the sequence.
        """
        return self.sequence
    
    def getLength(self) :
        """Determines the length of the sequence.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        int
            return a int corresponding to the length of the sequence.
        """
        return self.length

    def getTemperature(self) :
        """Determines the temperature of the conformation.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        int
            return a int corresponding to the temperature.
        """
        return self.temperature
    
    def getEnergy(self) :
        """Determines the energy of the conformation.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        int
            return a int corresponding to the energy.
        """
        return self.energy

    def setTemperature(self, temperature) :
        """Apply the given temperature.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        self.temperature = temperature

class MovementError(Exception):
    """Raised when a movement is impossible"""
    pass
