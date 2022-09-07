
class Residu :

    def __init__(self, residu=None, previous=None, next=None, coordinate=None) :
        self.residu = residu
        self.previous = previous
        self.next = next
        self.coordinate = coordinate

    def getResidu(self) :
        return self.residu

    def getPrevious(self) :
        return self.previous

    def getNext(self) :
        return self.next

    def getCoordinate(self) :
        return self.coordinate

    def setResidu(self, residu) :
        self.residu = residu

    def setPrevious(self, previous) :
        self.previous = previous

    def setNext(self, next) :
        self.next = next

    def setCoordinate(self, coordinate) :
        self.coordinate = coordinate

    def __str__(self) :
        return self.residu if self.residu else " "

