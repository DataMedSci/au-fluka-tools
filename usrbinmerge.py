import argparse
from glob import glob

import pandas
import numpy as np

def merge(inputfiles, outputfile):

    first_file = pandas.read_csv(inputfiles[0], sep=' ', names=('x','y','z','v'))
    result = first_file

    for file in inputfiles[1:]:
        data = pandas.read_csv(file, sep=' ', names=('x','y','z','v'))
        result['v'] += data['v']

    result['v'] /= len(inputfiles)

    result.to_csv(outputfile, sep=' ', header=False, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pattern', action='store', nargs='+',
                        help='input files (list of files or pattern')

    parser.add_argument('--outputfile', action='store', default='out_ascii.dat',
                        help='output file')

    args = parser.parse_args()


    inputfiles = args.pattern

    print("Input files: ", inputfiles)
    print("Output file: ", args.outputfile)


    merge(inputfiles=inputfiles, outputfile=args.outputfile)
