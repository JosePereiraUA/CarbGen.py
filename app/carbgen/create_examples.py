#!/usr/bin/env python

from residue import Residue
from sys import argv as args
import json

if __name__ == '__main__':
    """
    Example 1: VPC microcrystalite model
    Example 2: Single layer of functionalized oxide-graphene
    Example 3: Single non-functionalized graphene layer
    """

    examples = []
    x_list   = [6 , 10, 10]
    y_list   = [6 , 10, 10]
    z_list   = [3 , 1 , 1]
    pf_list  = [7 , 10, 1]
    co_list  = [9 , 9 , 0]
    coo_list = [3 , 5 , 0]
    ror_list = [18, 25, 0]
    ch_list  = [9 , 50, 0]

    for i in xrange(len(x_list)):
        residue = Residue(x_list[i], y_list[i], z_list[i], 'crv', 0.8 - 0.02 * float(pf_list[i]))
        residue.add_functionalizations(co_list[i], coo_list[i], ror_list[i], ch_list[i])
        residue.order_indexes()
        output = residue.printPDB_to_file('crv%d.pdb'%(i))
        molecule = {}
        molecule = {
            'x'  : x_list[i],
            'y'  : y_list[i],
            'z'  : z_list[i],
            'pf' : 0.8 - 0.02 * float(pf_list[i]), 
            'co' : co_list[i],
            'coo': coo_list[i],
            'ror': ror_list[i],
            'ch' : ch_list[i],
            'mol': "INSERT MOL HERE"
        }
        examples.append(molecule)

        #print ', '.join(json.dumps(molecule) for molecule in examples)
    print "%s%s%s" % ('[', ',\n'.join(json.dumps(molecule, sort_keys=True, indent=4, separators=(',', ': ')) for molecule in examples), ']')
