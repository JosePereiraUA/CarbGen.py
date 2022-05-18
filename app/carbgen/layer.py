from atom import Atom
from perlin import PerlinNoiseFactory
from functionalGroup import Storage
import numpy as np
import copy
import math

#import gc

class Layer:
    """
    Layer class encopasses all necessary variables and functions necessary to create and functionalize
    a layer of ATOMS;
    One RESIDUE is composed of one or more LAYERS.
    """
    def __init__(self):
        self.atoms = {}                                                    # Array of atoms (Dictionary format allows fast selection and is random by default)
        self.restoration_list = {'outer_carbons': [], 'inner_carbons': []} # List of atoms to restore after functional group addition


    """
    Pratical functions for iteration over the layer atoms
    """
    def __len__(self):
        count = 0
        for atom in self:
            count += 1
        return count


    def __iter__(self):
        for atom in self.atoms.iterkeys():
            if not atom.is_dead:
                yield atom


    def __getitem__(self, value):
        """
        Allows the correct picking of the carbon types inside the layer
        """
        if value == 'outer_carbons': return self.outer_carbons
        elif value == 'inner_carbons': return self.inner_carbons
        else: raise AttributeError("Layers has no {} attribute.".format(value))


    def createNewLayer(self, x, y, r = 1.4):
        """
        Creates a new 2 atom template and define the correct displacement vector used to copy the template;
        Outputs a fully filled atom layer.
        """
        self.x = x      # Number of carbon rings in x axis
        self.y = y      # Number of carbon rings in y axis
        
        #1. Define starting positions of the 2 atom template
        x1 = 0
        y1 = r
        x2 = r * math.cos(math.pi / 6.0)
        y2 = r * math.sin(math.pi / 6.0)

        #2. Define movement vectors to copy the 2 atom template
        ax, ay = 2 * x2, 0.0
        bx, by = x2, 1.5 * r

        #3. Populate layer with atoms by copying the original 2 atoms template
        step = 0
        atoms_list = []
        for j in range(self.y+1):
            connections_performed = 0
            for n, i in enumerate(range(-j // 2, self.x - j // 2)):
                dx = i * ax + j * bx
                dy = i * ay + j * by
                if (not self.y % 2 == 0) and (j == self.y) and (n == 0):
                    pass
                else:
                    atom1 = Atom(step, "C", [x1 + dx, y1 + dy, 0.0], 'ca')
                    step += 1
                atom2 = Atom(step, "C", [x2 + dx, y2 + dy, 0.0], 'ca')
                if j > 0 and connections_performed < self.x:
                    penalty = 1 if (j == self.y) and (not self.y % 2 == 0) else 0
                    if step == last_length + (1 + 2 * n) - penalty:
                        atom2.add_connection_to(atoms_list[step - (self.x * 2) - 1 + penalty])
                        connections_performed += 1
                step += 1
                if (not self.y % 2 == 0) and (j == self.y) and (n == 0):
                    pass
                else:
                    atom1.add_connection_to(atom2)
                    self.atoms[atom1] = None
                    atoms_list.append(atom1)
                if not n == 0:
                    atom1.add_connection_to(previous_atom)
                self.atoms[atom2] = None
                atoms_list.append(atom2)
                previous_atom = atom2
            last_length = step
        
        #4. Fill missing atoms on edges
        filled_atoms = 0
        for n, j in enumerate(range(0, self.y, 2)):
            penalty = 1 if j == self.y-1 else 0
            dx = (self.x - j / 2) * ax + j * bx
            dy = (self.x - j / 2) * ay + j * by
            atom1 = Atom(step, "C", [x1 + dx, y1 + dy, 0.0], 'ca')
            atom2 = (2 * self.x) + (2 * (self.x + 1)) * 2 * n - 1 - (n * 2) + filled_atoms
            atom3 = (2 * self.x) + (2 * (self.x + 1)) + (2 * (self.x + 1)) * (2*n) - 1 - penalty - (n * 2) + filled_atoms
            step += 1
            atom1.add_connection_to(atoms_list[atom2])
            atom1.add_connection_to(atoms_list[atom3])
            self.atoms[atom1] = None
            atoms_list.insert(atom2 + 1, atom1)
            filled_atoms += 1
        for n, j in enumerate(range(2, self.y + 1, 2)):
            dx = (-j / 2 - 1) * ax + j * bx
            dy = (-j / 2 - 1) * ay + j * by
            atom1 = Atom(step, "C", [x2 + dx, y2 + dy, 0.0], 'ca')
            atom2 = (2 * self.x) + (2 * (self.x + 1)) * 2 * n + 1
            atom3 = (2 * self.x) + (2 * (self.x + 1)) + (2 * (self.x + 1)) * (2*n) + 1
            step += 1
            atom1.add_connection_to(atoms_list[atom2])
            atom1.add_connection_to(atoms_list[atom3])
            self.atoms[atom1] = None
            atoms_list.insert(atom3, atom1)
        del atoms_list
        return self


    def createPores(self, pore_fraction, clean_sweeps = 15):
        """
        Create pores on the fully filled atom layer. Attaches a noise variable to each atom, based
        on the atom xy position. Noise is generated using Perlin Noise. A threshold (pore_fraction) is set;
        If the atom's noise value (between 0 and 1) is above the defined threshold,
        the atom is eficiently removed. Hanging atoms (with 1 connection) are removed in n clean sweeps.
        """

        #1. Create Perlin Noise
        perlin = PerlinNoiseFactory(2, octaves=3)

        #2. Create micropores based on Perlin Noise
        for atom in self:
            #2.1 Define noise value based on atom's xy position
            noise = perlin(atom.xyz[0], atom.xyz[1])
            #2.2 If the defined noise is above the query threshold, delete atom
            if noise > pore_fraction:
                atom.delete_connections()
                atom.is_dead = True

        #3. Clean rough edges and loose atoms
        for i in range(clean_sweeps):
            atoms_deleted = 0
            for atom in self:
                #Note: Loose atoms are defined as having only one connection
                if len(atom.connects) <= 1:
                    atom.delete_connections()
                    atom.is_dead = True
                    atoms_deleted += 1
            #Note: if no atoms were deleted in the previous sweep, break from process
            if atoms_deleted == 0: break
        return self


    def translate(self, v):
        """
        Move every atom of the current layer by an xyz vector.
        """
        for atom in self.atoms:
            atom.xyz = atom.xyz + v
        return self


    def renumber(self, start=0):
        """
        Renumber every atom of the current layer to start in a defined value.
        """
        for atom in self:
            atom.n = start
            start += 1 
        return self
    

    def set_carbon_types(self):
        """
        Divide the current layer's atoms based on connectivity:
        Outer carbons have 2 connections.
        Inner carbons have 3 connections.
        """

        #1. Initialize dictionary variables
        self.outer_carbons = {}
        self.inner_carbons = {}

        #2. Divide the layer's atoms based on connectivity
        for atom in self:
            if len(atom.connects) == 2:
                self.outer_carbons[atom] = None
            elif len(atom.connects) == 3:
                self.inner_carbons[atom] = None


    def restore_lists(self):
        """
        Return atoms lists to a previsouly saved state.
        """

        #1. Add all saved atoms in restorations lists to the layers inner and outer carbon lists
        for atom in self.restoration_list['outer_carbons']:
            self.outer_carbons[atom] = None
        for atom in self.restoration_list['inner_carbons']:
            self.inner_carbons[atom] = None

        #2. Reset saved state
        self.restoration_list['outer_carbons'] = {}
        self.restoration_list['inner_carbons'] = {}


    def add_ether_to(self, atom):
        """
        Turn a carbon atom into an oxygen atom, in an ether functional group.
        """

        #1. Place oxygen and remove it from avaliable list
        atom.symbol = 'O'
        self.outer_carbons.pop(atom)
        
        #2. Define correct charges to functional group and attached carbons
        excessCharge = Storage.FunctGroups['ror'].excessCharge
        split = excessCharge / len(atom.connects)
        
        #3. Add correct charge to functional oxygen
        atom.charge -= excessCharge
        
        #4. Remove neighbours from avaliable list and set correct charges to attached carbons
        restoration_list = self.restoration_list['outer_carbons']
        for connected_atom in atom.connects:
            connected_atom.charge += split
            if connected_atom in self.outer_carbons:
                self.outer_carbons.pop(connected_atom)
                #Note: Save state before functionalization
                restoration_list.append(connected_atom)


    def add_functionalization_to(self, key, atom, side_orientation, down_orientation, current_n):
        """
        Correctly orientate and place chosen functionalizations to the defined atom.
        """

        #1. Copy information regarding functional group
        fcn = copy.deepcopy(Storage.FunctGroups[key])

        #2. Set correct numbering based on last atom in layer
        fcn.renumber(current_n)

        #3. Set atomtypes and orientation
        atom.atmtype = 'c3'
        if side_orientation:
            atom.atmtype = 'ca'
            if key == 'co': atom.atmtype = 'c'
            fcn.sideways(atom, atom.connects.keys())
        elif down_orientation:
            fcn.downwards()

        #4. Place functionalization on correct atom
        fcn.move(np.asarray(atom.xyz))

        #5. Add connection to the chosen atom
        fcn[0].add_connection_to(atom)

        #6. Set charges
        atom.charge += fcn.carbonCharge
        split = fcn.excessCharge / len(atom.connects)
        for neighbour in atom.connects.keys(): neighbour.charge += split
        
        #7. Add placed atoms to list
        for fcn_atom in fcn.atoms:
            self.atoms[fcn_atom] = None
        
        #8. Remove current atom from avaliable list
        if side_orientation: self.outer_carbons.pop(atom)
        else: self.inner_carbons.pop(atom)
        
        #9. If the functionalization is a carboxyl, remove neighbouring atoms from avaliable list
        if key == 'coo':
            for neighbour in atom.connects.keys():
                if not side_orientation:
                    if neighbour in self.inner_carbons:
                        self.inner_carbons.pop(neighbour)
                        self.restoration_list['inner_carbons'][neighbour] = None