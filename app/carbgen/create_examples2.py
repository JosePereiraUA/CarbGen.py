#!/usr/bin/env python

from residue import Residue
from sys import argv as args
import json

if __name__ == '__main__':
    
    with open(args[1], 'r') as fin:
        data = "".join(fin.readlines())
        print "%s%s%s" % ('[', json.dumps(data), ']')