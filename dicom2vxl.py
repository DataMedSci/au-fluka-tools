# -*- coding: utf-8 -*-
from functools import reduce

from scipy.stats import itemfreq
import dicom, glob, struct
from numpy import *

from optparse import OptionParser


def ascii2vxl(filename, output, title, dims, sep=" ", hu_min=-1000, roworder='C'):
    data = fromfile(filename, sep=" ", dtype=short).reshape(dims, order=roworder)
    writect(output, title, data, hu_min, 1, 2, 1)


def numdicoms2vxl(basename, nmin, nmax, output, title, hu_min=-1000):
    """Convenience method for dicoms2vxl
       Takes a basename and 2 integers, and creates vxl file from the intermediate files.
       Example: num2dicoms("test",1,4,"ct.vxl","Some ct image") will create the ct.vxl file
       from test001.dcm, test002.dcm, test001.dcm and test003.dcm"""
    numbers = [str(x) for x in range(nmin, nmax + 1)]
    numbers = [("000")[0:(3 - len(x))] + x + ".dcm" for x in numbers]
    files = [basename + x for x in numbers]
    dicoms2vxl(files, output, title, hu_min)


def dicoms2vxl(files, output, title, hu_min=-1000):
    dicoms = map(lambda f: dicom.read_file(f), files)
    arrays = map(lambda ds: ds.pixel_array, dicoms)
    s = shape(arrays[0])
    array = reduce(lambda x, y: append(x, y), arrays)
    array = array.reshape(s[0], s[1], len(arrays))
    ds = dicoms[0]
    dx, dy = ds.PixelSpacing
    dz = ds.SliceThickness
    writect(output, title, array, -1000, dx, dy, dz)


def dicom2vxl(filename, output, title, hu_min):
    ds = dicom.read_file(filename)
    array = ds.pixel_array
    dx, dy = ds.PixelSpacing
    dz = ds.SliceThickness

    writect(output, title, array, -1000, dx, dy, dz)


def writect(filename, title, a, hu_min, dx, dy, dz):
    f = open(filename, "wb")
    ct = a - hu_min  # Look Ma! This is about 10 lines of Fortran
    mo = amax(ct)
    ireq = set()
    kreq = zeros(mo, dtype=int16)
    print("Making data for VXL file")

    no = 0
    for x in ct.flat:
        if not (x in ireq):
            kreq[x - 1] = no
            ireq.add(x)
            print("HU:" + str(x) + " FLUKA val:" + str(no))
            no += 1

    dims = shape(a)
    # Danger! Writing Fortran data structures. Do NOT try at home!
    print("Writing FLUKA-Fortran data")
    for x in range(80 - len(title)):
        title += " "

    f.write(struct.pack("=i80si", 80, title, 80))  # WRITE TITLE
    f.write(struct.pack("=iiiiiii", 20, dims[0], dims[1], dims[2], no - 1, mo, 20))
    f.write(struct.pack("=idddi", 24, dx, dy, dz, 24))
    f.write(struct.pack("=i", size(ct) * 2))
    ctshape = shape(ct)
    # Writing 3d array Fortran style.
    for x in range(ctshape[0]):
        for y in range(ctshape[1]):
            for z in range(ctshape[2]):
                f.write(struct.pack("=h", ct[x, y, z]))
    f.write(struct.pack("=i", size(ct) * 2))
    f.write(struct.pack("=i", size(kreq) * 2))
    kreq.tofile(f)
    f.write(struct.pack("=i", size(kreq) * 2))


if __name__ == "__main__":
    version = "0.1"

    parser = OptionParser()
    parser.add_option("-t", "--title", dest="title",
                      default="",
                      help="FLUKA title of the voxel geometry", metavar="string")

    parser.add_option("--hu_min", dest="hu_min",
                      default=-1000,
                      help="Minimum HU unit. Default: -1000")
    parser.add_option("-o", "--output", dest="output",
                      help="Name of the output file. Example ct.vxl ")
    (options, args) = parser.parse_args()

    filename = args[0]
    filenames = glob.glob(filename)
    if len(filenames) == 1:
        dicom2vxl(filenames[0], options.output, options.title, options.hu_min)
        print("Vxl file written")
    elif len(filenames) == 0:
        print("Please enter more than 1 filename")
    else:
        dicoms2vxl(filenames, options.output, options.title, options.hu_min)
        print("Vxl file written")
