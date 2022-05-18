from StringIO import StringIO
import zipfile
import os

from residue import Residue


def validate_parameters(*args, **kwargs):
    return True


def batch_create(n_batches, prefix, x, y, z, pf, co, coo, ror, ch, append_dir=None):
    """
    Creates a batch of n residues, with the given input parameters;
    Adds necessary files to an output zip file.
    """
    
    #1. Create memory buffer
    mem_buffer = StringIO()
    
    #2. Create zip file for output
    zipf = zipfile.ZipFile(mem_buffer, 'w', zipfile.ZIP_DEFLATED)
    
    #3. Add forcefield folder to zip file
    if append_dir is not None and os.path.isdir(append_dir):
        for root, dirs, files in os.walk(append_dir):
            for file in files:
                arcname = os.path.join(
                                    os.path.basename(append_dir),
                                    os.path.relpath(root, append_dir),
                                    file)
                zipf.write(os.path.join(root, file), arcname=arcname)
                

    #4. Create and add residues to zip file
    for n in xrange(n_batches):
        #4.1 Create residue
        residue = Residue(x, y, z, prefix, pf)
        #4.2 Add functionatilzations to residue
        residue.add_functionalizations(co, coo, ror, ch)
        #4.3 Correctly order atom indexes in residue
        residue.order_indexes()
        #4.4 Print requested file formats to memory buffer
        for ftype in ('MOL2', 'PDB','TOP','INF'):
            f = getattr(residue, 'print' + ftype)
            buffer = StringIO()
            f(buffer)
            fname = '{prefix}_{n}.{extension}'.format(prefix=prefix, n=n, extension=ftype.lower())
            #4.5 Add output files to zip file
            zipf.writestr(fname, buffer.getvalue())
    zipf.close()
    
    mem_buffer.seek(0)
    content = mem_buffer.read()
    mem_buffer.close()

    return content



if __name__ == '__main__':
    from time import time as stopwatch
    import argparse

    """
    Offline CarbGen Tool:
    This python script is intended to help the test and optimization of the CarbGen tool.
    """

    start = stopwatch()

    #1. Parse arguments
    main_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main_parser.add_argument('-x', '--x_size', default=6, type=int, action='store',
                    help='Number of carbon rings in x axis')
    main_parser.add_argument('-y', '--y_size', default=6, type=int, action='store',
                    help='Number of carbon rings in y axis')
    main_parser.add_argument('-z', '--z_size', default=2, type=int, action='store',
                    help='Number of carbon layers')
    main_parser.add_argument('-co', '--carbonyl', default=9.0, type=float, action='store',
                    help='Carbonyl percentage')
    main_parser.add_argument('-coo', '--carboxyl', default=3.0, type=float, action='store',
                    help='Carboxyl percentage')
    main_parser.add_argument('-ror', '--ether', default=18.0, type=float, action='store',
                    help='Ether percentage')
    main_parser.add_argument('-ch', '--hydro', default=9.0, type=float, action='store',
                    help='Hydrogenation percentage')
    main_parser.add_argument('-n', '--name', default="crv", type=str, action='store',
                    help='Residue name')
    main_parser.add_argument('-pf', '--pore_fraction', default=0.6, type=float, action='store',
                    help='Pore fraction')
    main_parser.add_argument('-b', '--batch', default=1, type=int, action='store',
                    help='Number of batch processes')
    args = main_parser.parse_args()

    #2. Create residues in batch
    residues = batch_create(args.batch, args.name, args.x_size, args.y_size, args.z_size, args.pore_fraction,
        args.carbonyl, args.carboxyl, args.ether, args.hydro)

    #3. Save zip file to output zip folder
    with open('{name}.zip'.format(name = args.name), 'wb') as fout:
        fout.write(residues)

    end = stopwatch()
    print "\n\033[94mProduced residue in {time:.2f}s.\033[0m".format(time = end - start)

    