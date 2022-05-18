class Atom:
    """
    Atom class defines a point in the LAYER, with a defined position, charge, atomtype and connections
    """
    def __init__(self, n="NaN", symbol='NaN', xyz=[0.,0.,0.], atmtype='NaN', charge = 0.000):
        self.n = n                      # Atom number 
        self.symbol = symbol            # Atom symbol
        self.xyz = xyz                  # Atom position
        self.connects = {}              # Array of atoms connected to this atom
        self.atmtype = atmtype          # Atom type
        self.charge = charge            # Atom charge
        self.is_dead = False            # An atom defined as "dead" is not accounted for in iterations


    def move(self, v, n):
        """
        Copies the current atom, to a new position based on an xyz vector, and sets correct atom number.
        """
        new = ([self.xyz[0]+v[0], self.xyz[1]+v[1], self.xyz[2]+v[2]])
        return Atom(self.n+n, self.symbol, new, self.atmtype, self.charge)


    def add_connection_to(self, atom):
        """
        Adds this atom to the target connection list and vice-versa.
        """
        self.connects[atom] = None
        atom.connects[self] = None


    def delete_connections(self):
        """
        Delete the record of this atoms in all neighbouring atoms connect list.
        """
        for atom in self.connects:
            atom.connects.pop(self, None)
