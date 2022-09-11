class Residu :

    def __init__(self, residu=None, previous=None, next=None, coordinate=None, hydrophobicity=None) :
        self.residu = residu
        self.previous = previous
        self.next = next
        self._coordinate = coordinate
        self.hydrophobicity = hydrophobicity

    def translateHP(self) :
        hydrophobicity = {"A": "H", "I": "H", "L": "H", "M": "H", "F": "H", 
        "W": "H", "Y": "H", "V": "H", "G": "H", "P": "H",
        "T": "P", "K": "P", "R": "P", "H": "P", "D": "P", "E": "P", "S": "P", 
        "N": "P", "Q": "P", "C": "P", "U": "P"}
        if (self.residu in hydrophobicity) :
            self.hydrophobicity = hydrophobicity[self.residu]

    def getResidu(self) :
        return self.residu

    def getHydrophobicity(self) :
        return self.hydrophobicity

    def getPrevious(self) :
        return self.previous

    def getNext(self) :
        return self.next

    def getCoordinate(self) :
        return self.coordinate

    def setResidu(self, residu) :
        self.residu = residu
    
    def setHydrophobicity(self, hydrophobicity) :
        self.hydrophobicity = hydrophobicity

    def setPrevious(self, previous) :
        self.previous = previous

    def setNext(self, next) :
        self.next = next

    def setCoordinate(self, coordinate) :
        self.coordinate = coordinate

    def __str__(self) :
        return self.residu if self.residu else " "

