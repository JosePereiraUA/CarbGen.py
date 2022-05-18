from conversion import residue_to_gmx
from functionalGroup import Storage
from layer import Layer
import numpy as np
import random

class Residue: 
    """
    Residue class encopasses all necessary functions to create and functionalize a new residue
    """
    def __init__(self, x, y, z, name, pore_fraction):
        """
        A residue is compromised of one or more LAYERS, each made of ATOMS;
        By default, creates a new blank residue with given size, number of layers and pore fraction
        """
        self.x = x                      # Number of carbon rings in x axis
        self.y = y                      # Number carbon rings in y axis
        self.z = z                      # Number of carbon layers
        self.name = str(name)           # Residue name, used in printMOL2 function
        self.functional_groups = {'co': 0, 'coo': 0, 'ch': 0, 'ror': 0}

        self.createNewResidue(x, y, z, pore_fraction = pore_fraction)
    

    """
    Pratical functions for iteration over the different layers of a residue
    """
    def __len__(self):
        count = 0
        for atom in self:
                count += 1
        return count


    def __iter__(self):
        for layer in self.layers:
            for atom in layer:
                yield atom


    def order_indexes(self):
        """
        ATOMS in LAYERS are organized in dictionaries (no order);
        Creates a new variable self.atoms that lists all atoms in the residue, ordering by atom INDEX;
        Used in order to print output files
        """
        self.atoms = [atom for atom in self]
        self.atoms.sort(key=lambda atom: atom.n)


    def createNewResidue(self, x, y, z, pore_fraction = 1, r = 1.4, spacing = 3.4):
        """
        Creates a new black residue, with a given size (x * y) and a number of layers (z)
        Adds pores, based on the given pore_fraction
        """
        
        #1. Initialize variables
        self.layers = []                                         #Array of layers in the residue
        displ_vec = np.array([r*np.cos(np.pi/6), r/2, spacing])  #Displacement vector used to translate new layers
        atoms_per_layer = 0                                      #Atoms per layer get increased in order to keep correct numbering of atoms INDEX since, because of pores, layers can have different number of atoms
        
        #2. Create necessary layers
        for i in range(z):
            
            #2.1 Create new layer, and translate in x, y and z for correct pi-pi stacking
            new_layer = Layer()
            new_layer.createNewLayer(x, y).translate(displ_vec*i)
            
            #2.2 Add pores and renumber the layer to start on the last index of previous layers
            new_layer.createPores(pore_fraction).renumber(atoms_per_layer)
            
            #2.3 Detect if the carbons are in edges (2 connects) or centers (3 connects)
            new_layer.set_carbon_types()
            
            #2.4 Add the last index of the layer atoms_per_layer
            atoms_per_layer += len(new_layer)
            
            #2.5 Append layer to residue list
            self.layers.append(new_layer)


    def add_functionalizations(self, co = 0, coo = 0, ror = 0, ch = 0):
        """
        Add functionalizations IN ORDER to the blank residue;
        CO  = Carbonyls;
        COO = Carboxyls;
        ROR = Ethers;
        CH  = Hydrogen terminations
        """
        #1. initialize variables for each of the functionalizations currently supported
        self.ror = ror
        self.co = co
        self.coo = coo
        self.ch = ch

        #2. Exit function if no functionalizations were requested
        if (co == 0) and (coo == 0) and (ror == 0) and (ch == 0):
            return self

        #3. Add functionalizations in order ROR -> COO -> CO -> CH
        self.add_functionalization('ror', ror)
        self.restore_lists()                    #  During the functionalization process, certain positions
        self.add_functionalization('coo', coo)  # are removed from the possible anchor points, given
        self.restore_lists                      # chemical constrains. However, those constrains do not
        self.add_functionalization('co', co)    # exist for other functional group types. Possibility lists
        self.add_functionalization('ch', ch)    # are therefore restored before new functionalization.


    def restore_lists(self):
        """
        Restore possible anchor points lists to defaults
        """
        for layer in self.layers:
            layer.restore_lists()


    def add_functionalization(self, key, value):
        """
        Place SPECIFIC functionalization in the residue.
        Control amount placed to achieve desired quantity.
        Choose random atom to place.
        Orient functionalization correctly.
        """
        
        #1. Exit function if no functionalization of this type was requested.
        if value == 0: return
        
        #2. Count total carbon atoms to control correct quantification of functional groups attached
        totalC = self.count('C')
        
        #3. Count existing atoms for correct numbering of newly attached functional groups
        total = len(self)
        
        #4. Initialize variables to control correct quantification of functional groups
        percentage = 0
        fcn_count = 0
        
        #5. Keep adding functionalizations until placed percentage = requested percentage
        while percentage < value:
            
            #5.1 Initialize orientation control booleans
            side_orientation = False
            down_orientation = False #Note: Default is up orientation
            
            #5.2 Choose carbon types that this functionalization can be attached to
            layers_to_search = 'outer_carbons' # All functionalizations can be placed on outer_carbons
            if key == 'coo':                   # COO functionalizations can aditionally be attached to inner_carbons
                layers_to_search = random.choice(['outer_carbons', 'inner_carbons'])
                try:
                    #Possible layers include all layers with exisitng free eligible spots for functionalization 
                    possibleLayers = []
                    #COO CASE: Search outer_ and inner_carbons for first and last layer
                    for i in [0, -1]:
                        if len(self.layers[i][layers_to_search].keys()) > 0:
                            possibleLayers.append(self.layers[i])
                    #COO CASE: If existent, search middle layers for outer_carbons ONLY
                    if len(self.layers) > 2:
                        for i in xrange(1, len(self.layers) - 1):
                            if len(self.layers[i]['outer_carbons'].keys()) > 0:
                                possibleLayers.append(self.layers[i])
                    
                    #5.3 Pick random layer to attach functionalization
                    random_layer = random.choice(possibleLayers)
                    del possibleLayers
                except IndexError:
                    #5.4 If all eligible spots are exhausted, break the 'while' cycle
                    break
            else:
                try:
                    #GENERAL CASE: Search all layers for outer_carbons ONLY
                    #5.3 Pick random layer to attach functionalization
                    random_layer = random.choice([layer for layer in self.layers if len(layer[layers_to_search].keys()) > 0])
                except IndexError:
                    #5.4 If all eligible spots are exhausted, break the 'while' cycle
                    break
            
            #5.6 Pick random atom from the picked random layer
            #In the case of COO functionalization in middle layers, only search outer_carbons
            if (len(self.layers) > 2) and key == 'coo' and random_layer in self.layers[1:(len(self.layers) - 1)]:
                    random_atom = random.choice(random_layer['outer_carbons'].keys())
            #In the general case, search both inner_ and outer_carbons
            else:
                random_atom = random.choice(random_layer[layers_to_search].keys())
            
            #5.7 Correctly orientate functionalization (Note: Default is up orientation)
            if key == 'ch' or key == 'co':
                side_orientation = True #CH and CO functionalizations are ALWAYS side-oriented
            elif random_atom in random_layer.outer_carbons:
                side_orientation = True #COO functionalizations in outer_carbons are side-oriented
            elif len(self.layers) == 1 and random.choice([True, False]):
                down_orientation = True #If only 1 layers exists, randomly pick top or bottom orientation
            elif len(self.layers) > 1 and random_atom in self.layers[0]:
                down_orientation = True #If multiple layers exists and functionalization is attached to
                                        #first layer, use down-orientation
            
            #5.8 Place functionalization
            if key == 'ror':
                random_layer.add_ether_to(random_atom)
                totalC -= 1
            else:
                random_layer.add_functionalization_to(key, random_atom, side_orientation, down_orientation, total)
            
            #5.9 Count number of carbons involved in attachment for functional content calculations in INFO file
            fcn_count = Storage.FunctGroups[key].carbon_atoms_attached
            
            #5.10 Add the number of atoms added to total count, for correct numbering of future functionalizations
            total += len(Storage.FunctGroups[key].atoms)
            
            #5.11 Add the number of atoms added to self.functional_groups for functional content calculations in INFO file
            self.functional_groups[key] += fcn_count
            
            #5.12 Calculate placed percentage
            percentage = (self.functional_groups[key] / totalC) * 100


    def count(self, atom_element = None):
        """
        Count atoms. If an element is specified, count only atoms of the query element.
        """
        count = 0
        for atom in self:
            if atom_element == None: count += 1
            elif atom.symbol == atom_element: count += 1
        return float(count)


    def printPDB(self, filout):
        """
        Print residue information to PDB file format. filout variable should be a StringIO buffer.
        """
        
        #1. Write ATOM records
        for atom in self.atoms:
                form = 'ATOM %6d  %-3s UNK   %3d    %8.3f%8.3f%8.3f ; %-3s %-3.4f\n'
                filout.write(form % (atom.n, atom.symbol, 1, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.atmtype, atom.charge))
        
        #2. Write CONECT records
        filout.write('TER\n')
        for atom in self.atoms:
                filout.write("CONECT%5d" % atom.n)
                filout.write(''.join('%5d' % (j.n, ) for j in atom.connects))
                filout.write("\n")

    
    def printPDB_to_file(self, file_name):
        """
        Print residue information to PDB file format. filout variable should be a file name.
        """
        
        with open(file_name, 'w') as filout:
            #1. Write ATOM records
            for atom in self.atoms:
                form = 'ATOM %6d  %-3s UNK   %3d    %8.3f%8.3f%8.3f ; %-3s %-3.4f\n'
                filout.write(form % (atom.n, atom.symbol, 1, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.atmtype, atom.charge))
            
            #2. Write CONECT records
            filout.write('TER\n')
            for atom in self.atoms:
                filout.write("CONECT%5d" % atom.n)
                filout.write(''.join('%5d' % (j.n, ) for j in atom.connects))
                filout.write("\n")


    def printMOL2(self, filout):
        """
        Print residue information in MOL2 file format. filout variable should be a StringIO buffer.
        """
        
        #1. Name residue
        resN = "R" + self.name
        
        #2. Define correct bonds in MOL2 format
        bonds = self.mol2Bonds()
        
        #3. Print HEADER
        filout.write("@<TRIPOS>MOLECULE\n%s\n%d %d %d\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" % (resN, len(self), len(bonds), 1))
        
        #4. Print ATOM records
        form  = '%-6d %-4s %6.3f %6.3f %6.3f %-9s %2d %-8s %9.6f\n'
        for atom in self.atoms:
            filout.write(form % (atom.n + 1, atom.symbol, atom.xyz[0], atom.xyz[1], atom.xyz[2], atom.atmtype, 1, 'CRV1', atom.charge))
        
        #5. Print BOND records
        filout.write("@<TRIPOS>BOND\n")
        form = '%d %d %d %d\n'
        for n, bond in enumerate(bonds):
            filout.write(form % (n, bond[0], bond[1], 1))
        
        #6. Print FOOTER
        filout.write("@<TRIPOS>SUBSTRUCTURE\n%-3d %-7s %-3d %-5s %3d %-5s %s\n" % (1, resN, 1, 'GROUP', 1, '****', resN))


    def mol2Bonds(self):
        """
        Creates a list of atom pairs connected, with correct numbering for MOL2 format.
        """
        bonds = []
        for atom1 in self.atoms:
            for atom2 in atom1.connects:
                #Note: To eliminate repetitions, only add bonds where the second atom number is greater
                if atom2.n > atom1.n: bonds.append([atom1.n + 1, atom2.n + 1])
        return bonds


    def printTOP(self, buffer):
        """
        Use conversion.py function "residue_to_gmx" to write residue information in topology format.
        """
        residue_to_gmx(self, buffer)


    def printINF(self, filout):
        """
        Calculate and print additional information regarding the residue chemistry to INFO file.
        """
        
        #1. Count the total number of atoms
        tot = len(self)
        
        #2. Perform elemental analysis
        elem = self.elementalAnalysis()
        
        #3. Complete functional content information with number of sp2 carbons
        self.functional_groups['sp2'] = self.count('C') - self.functional_groups['co'] - self.functional_groups['coo'] - self.functional_groups['ch'] - self.functional_groups['ror']
        
        #4. Define conversion table for correct elemental analysis in (w/w)
        w, totalW = {'C': 12.0107, 'O': 15.9994, 'H': 1.00794}, 0
        
        #5. Calculate residue molar mass
        for key in elem: totalW += elem[key]*w[key]
        
        #6. Print RESIDUE DETAILS records
        filout.write("[ Residue details ]\n")
        filout.write("Atom count: %d\n" % (tot))
        filout.write("Residue charge: %d\n" % (self.functional_groups['coo'] * -1))
        filout.write("Residue molar mass: %d g/mol\n" % (totalW))
        
        #7. Print ELEMENTAL ANALYSIS RECORDS
        filout.write("\n[ Elemental Analysis (mol)  ]\n")
        for key in elem:
            per = (float(elem[key])/tot)*100
            filout.write("%s : %3.2f%% ( %3d / %3d )\n" % (key, per, elem[key], tot))
        filout.write("\n[ Elemental Analysis (mass) ]\n")
        for key in elem:
            per = ((float(elem[key])*w[key])/totalW)*100
            filout.write("%s : %3.2f%%\n" % (key, per))
        
        #8. Print FUNCTIONAL CONTENT records
        filout.write("\n[ Functional content (mol)  ]\n")
        tot = self.count('C')
        for key in self.functional_groups:
            per = (float(self.functional_groups[key])/tot)*100
            if not self.functional_groups[key] == 0:
                filout.write("%-4s : %3.2f%% ( %3d / %3d )\n" % (key, per, self.functional_groups[key], tot))
        filout.write("\n---\n")


    def elementalAnalysis(self):
        """
        Count number of atoms of each element considered.
        """
        out = {'C': 0, 'O': 0, 'H': 0}
        for atom in self:
            out[atom.symbol] += 1
        return out


    def printJSON(self):
        """
        Prints residue information in JSON format,
        as described in https://web.chemdoodle.com/docs/chemdoodle-json-format/
        """
        #1. Initialize variables
        output = []
        molecule = {'a': [], 'b': []}
        
        #2. Correctly define bonds
        bonds = self.mol2Bonds()
        
        #3. Add ATOM records
        for atom in self:
            new_atom = {'l': atom.symbol, 'x': atom.xyz[0], 'y': atom.xyz[1], 'z': atom.xyz [2]}
            molecule['a'].append(new_atom)
        
        #4. Add BOND records
        for bond in bonds:
            new_bond = {'e': bond[0] - 1, 'b': bond[1] - 1}
            molecule['b'].append(new_bond)
        
        output.append(molecule)
        return output

