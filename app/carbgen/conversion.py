import cStringIO

def gen_path(pivot, max_length, container, path=None):
    '''
    gen_path: utility function to generate all
    paths within a max_length of a pivot node
    '''
    if path is None:
        path = [pivot]
    
    if len(path) == max_length:
        path = tuple(path)
        if (path not in container) and (path[0].n < path[-1].n):
            container[path] = True
    
    else:
        for other in pivot.connects:
            if other not in path:
                gen_path(other, max_length, container, path + [other])
    

def residue_to_gmx(residue, buffer):

    atoms = [atom for atom in residue]
    atoms.sort(key=lambda atom: atom.n)
    for atom in atoms:
        atom.n += 1

    def node_list_to_str(nodes, funct):
        items = [atom.n for atom in nodes] + [funct]
        return '{0}  ;  {1}\n'.format(
            ''.join('%6d'%i for i in items),
            ' - '.join('%-5s'%(atom.symbol) for atom in nodes)
        )
        

    def write_entries(buffer, items, funct, *headers):
        entries = sorted(items, key=lambda item: [atom.n for atom in item])
        
        for header in headers:
            buffer.write(header + '\n')

        for entry in entries:
            s = node_list_to_str(entry, funct)
            buffer.write(s)
        return entries

    # -----------------------------------------------------
    # HEADER
    # -----------------------------------------------------
    # buffer = cStringIO.StringIO()
    buffer.write('; Include forcefield parameters\n')
    buffer.write('#include "sheet_ga.ff/forcefield.itp"\n')
    
    buffer.write('\n[ moleculetype ]\n')
    buffer.write('; name      nrexcl\n')
    buffer.write(' crv      3\n')
    
    # -----------------------------------------------------
    # ATOMS
    # -----------------------------------------------------
    buffer.write('\n[ atoms ]\n')
    buffer.write(';     nr     type    resnr  residue     atom     cgnr   charge     mass\n')
    qtot = 0.0
    for atom in atoms:
        qtot += atom.charge
        s = '%8d %8s %8d %8s %8s %8d %8.4f %8.3f ; qtot %8.4f\n'%(atom.n, atom.atmtype, 1, 'crv', atom.symbol, atom.n, atom.charge, 1.0, qtot)
        buffer.write(s)


    # -----------------------------------------------------
    # BONDS
    # -----------------------------------------------------
    bonds = [(atom, other) for atom in residue for other in atom.connects if other.n > atom.n]
    write_entries(buffer, bonds, 1, '',
                    '[ bonds ]',
                    ';   ai    aj funct')

    # -----------------------------------------------------
    # ANGLES
    # -----------------------------------------------------
    container = {}
    for atom in atoms:
        gen_path(atom, 3, container)
    write_entries(buffer, container.keys(), 1, '',
                    '[ angles ]',
                    ';   ai    aj    ak funct')

    # -----------------------------------------------------
    # PROPER DIHEDRALS
    # -----------------------------------------------------
    container = {}
    for atom in atoms:
        gen_path(atom, 4, container)
    proper = write_entries(buffer, container.keys(), 3, '',
                    '[ dihedrals ] ; proper',
                    ';   ai    aj    ak    al funct')

    pairs = [(p[0], p[3]) for p in proper]
    write_entries(buffer, pairs, 1, '',
                    '[ pairs ]',
                    ';   ai    aj funct')

    # -----------------------------------------------------
    # IMPROPER DIHEDRALS
    # -----------------------------------------------------
    improper = []
    for atom in atoms:
        if len(atom.connects) == 3:
            a1,a2,a4 = sorted(atom.connects, key=lambda n: (n.atmtype, n.n))
            improper.append((a1, a2, atom, a4))
    write_entries(buffer, improper, 1, '',
                    '[ dihedrals ] ; improper',
                    ';   ai    aj    ak    al funct')
    
    # -----------------------------------------------------
    # SYSTEM & MOLECULES
    # -----------------------------------------------------
    buffer.write('\n[ system ]\n')
    buffer.write('crv\n')
    
    buffer.write('\n[ molecules ]\n')
    buffer.write('; Compound        nmols\n')
    buffer.write('crv      1\n')