"""Python file grouping the Residu class as well as the associated functions.
"""

class Residu :
    """A class representing a residu of a protein.

    Attributes
    ----------
        residu : string
            a string corresponding to a residu of a letter
        hydrophobic : string
            a string corresponding to a polarity
    """
    def __init__(self, residu=None, hydrophobic=None) :
        """The constructor of the Residu object

        Parameters
        ----------
        self : object
            a parameter in instance method.
        residu : string
            a parameter corresponding to a residu
        hydrophobic : string
            a parameter corresponding to a hydrophobic

        Returns
        -------
        None
            return nothing.
        """
        self.residu = residu
        self.hydrophobic = hydrophobic

    def translateHP(self) :
        """Translates the residue according to the HP model (hydrophobic or polar).

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        None
            return nothing.
        """
        hydrophobic = {"A": "H", "I": "H", "L": "H", "M": "H", "F": "H", 
        "W": "H", "Y": "H", "V": "H", "G": "H", "P": "H",
        "T": "P", "K": "P", "R": "P", "H": "P", "D": "P", "E": "P", "S": "P", 
        "N": "P", "Q": "P", "C": "P", "U": "P"}
        if (self.residu in hydrophobic) :
            self.hydrophobic = hydrophobic[self.residu]

    def getResidu(self) :
        """Determines the string of the residu.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        string
            return a string corresponding to the residu.
        """
        return self.residu

    def getHydrophobic(self) :
        """Determines the string of the polarity.

        Parameters
        ----------
        self : object
            a parameter in instance method.

        Returns
        -------
        string
            return a string corresponding to the polarity.
        """
        return self.hydrophobic

    def setResidu(self, residu) :
        """Apply the given residu.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        residu : string
            a parameter corresponding to a residu

        Returns
        -------
        None
            return nothing.
        """
        self.residu = residu
    
    def setHydrophobic(self, hydrophobic) :
        """Apply the given hydrophobic.

        Parameters
        ----------
        self : object
            a parameter in instance method.
        hydrophobic : string
            a parameter corresponding to a polarity

        Returns
        -------
        None
            return nothing.
        """
        self.hydrophobic = hydrophobic

    def __str__(self) :
        return self.residu if self.residu else " "

