from atom import Atom
import numpy as np

class FunctionalGroup:
    """
    Defines necessary functions and variables for handling a functional group.
    """
    def __init__(self, atoms=None, carbonCharge = 0, excessCharge = 0, carbon_atoms_attached = 0):
        if atoms==None: self.atoms = []
        else: self.atoms = atoms            # List of atoms involved in the functional group
        self.carbonCharge = carbonCharge    # Charge of the target carbon atom for attachment
        self.excessCharge = excessCharge    # Excess charge to be spread among neighbouring carbon atoms
        self.carbon_atoms_attached = carbon_atoms_attached


    def __getitem__(self, key):
        return self.atoms[key]


    def setAtom(self, atom):
        """
        Add an atom to the funcitonal group.
        """
        self.atoms.append(atom)


    def renumber(self, n):
        """
        Set current number based on an input int.
        """
        for atom in self.atoms:
            atom.n = n
            n += 1


    def downwards(self):
        """
        Move necessary atoms to invert the facing direction of the functional group.
        """
        for atom in self.atoms: atom.xyz[2] = - atom.xyz[2]


    def sideways(self, at0, neighbours):
        """
        Define neighbour atoms, calculate the rotation vector to orient the funcitonal group sideways.
        """

        #1. Define neighbouring atoms
        at1 = neighbours[0]
        at2 = neighbours[1]

        #2. Define rotation vector based on neighbouring atoms positions
        v1 = np.asarray(at1.xyz) - np.asarray(at0.xyz)
        v2 = np.asarray(at2.xyz) - np.asarray(at0.xyz)
        v3 = v1 + v2
        w = np.cross([0,0,1], v3)

        #3. Apply rotation vector to all atoms of the current functional group
        for fcnatom in self.atoms:
            new_xyz = self.axisRotation(w, -(np.pi / 2), np.asarray(fcnatom.xyz))
            fcnatom.xyz = new_xyz


    def axisRotation(self, w, t, v):
        """
        Rotate an xyz vector by a given rotation vector around a given axis.
        """
        m = 0
        try:
            n, m = v.shape
        except:
            n, = v.shape
        wn = w / np.sqrt(np.dot(w, w))
        ct = np.cos(t)
        if m:
            return v * ct + np.cross(wn, v) * np.sin(t) + wn * (1 - ct) * np.reshape(np.repeat((wn * v).sum(1), 3), (-1, 3))
        return v * ct + np.cross(wn, v) * np.sin(t) + wn * (1 - ct) * (wn * v).sum()


    def move(self, v):
        """
        Move all the atoms in the fucntional group by an xyz vector.
        """
        for atom in self.atoms:
            atom.xyz = atom.xyz + v
        return self


class Storage:
    """
    Stores information of functional groups regarding number of atoms, atom positions and charges.
    """
    FunctGroups = {}
    
    #1. Hydrogen terminations
    FunctGroups['ch'] = FunctionalGroup(carbonCharge = 0.16, carbon_atoms_attached = 1)
    FunctGroups['ch'].setAtom(Atom('NaN', 'H', [0.0, 0.0, 0.8], 'ha', -0.16))
    
    #2. Ethers
    FunctGroups['ror'] = FunctionalGroup(excessCharge = 0.25, carbon_atoms_attached = 2)
    
    #3. Carboxyls
    FunctGroups['coo'] = FunctionalGroup(carbonCharge = 0.03, excessCharge = -0.34, carbon_atoms_attached = 1)
    FunctGroups['coo'].setAtom(Atom('NaN', 'C',[ 0.00,  0.00, 1.49],'c2', 0.83))
    FunctGroups['coo'].setAtom(Atom('NaN', 'O',[ 1.04,  0.17, 2.10],'o', -0.76))
    FunctGroups['coo'][1].add_connection_to(FunctGroups['coo'][0])
    FunctGroups['coo'].setAtom(Atom('NaN', 'O',[-1.15, -0.19, 2.16],'o', -0.76))
    FunctGroups['coo'][2].add_connection_to(FunctGroups['coo'][0])
    
    #4. Carbonyls
    FunctGroups['co'] = FunctionalGroup(carbonCharge = 0.54, excessCharge = 0.02, carbon_atoms_attached = 1)
    FunctGroups['co'].setAtom(Atom('NaN', 'O', [0.0, 0.0, 1.03], 'o', -0.56))