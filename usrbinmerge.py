import argparse
from glob import glob

import pandas
import numpy as np


def merge(input_files, output_file):
    first_file = pandas.read_csv(input_files[0], sep=' ', names=('x', 'y', 'z', 'v'))
    result = first_file

    for file in input_files[1:]:
        data = pandas.read_csv(file, sep=' ', names=('x', 'y', 'z', 'v'))
        print(np.max(data['x']), np.max(data['y']), np.max(data['z']))
        result['v'] += data['v']

    result['v'] /= len(input_files)

    for axis in ('x', 'y', 'z'):
        if np.unique(result[axis]).size == 1:
            result.drop(axis, inplace=True, axis=1)

    result.to_csv(output_file, sep=' ', header=False, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pattern', action='store', nargs='+',
                        help='input files (list of files or pattern')

    parser.add_argument('--outputfile', action='store', default='out_ascii.dat',
                        help='output file')

    args = parser.parse_args()

    print("Input files: ", args.pattern)
    print("Output file: ", args.outputfile)

    merge(input_files=args.pattern, output_file=args.outputfile)
