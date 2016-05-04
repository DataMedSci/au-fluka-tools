#!/usr/bin/env python
#########################################################################
# Copyright (C) 2010 Niels Bassler, bassler@phys.au.dk
#
# This program is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation; either version 3 of the License, or (at your 
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along 
# with this program; if not, see <http://www.gnu.org/licenses/>.
#
# usrbin2ascii.py
#
# Converting usrbin files to ascii. 
#
#
# NB: Niels Bassler, bassler@phys.au.dk,
#
# TODO:
# - check if input files are binary before processing, improve file checking
# - allow usrbin2ascii foo*   instead of     usrbin2ascii "foo*"
# - less verbose option ?
# - merge usrbin2ascii wtih usrtrack2ascii
# - convert to shell scripts which sets flair path for python instead of this solution ?
# - multiple detectors

import os
import sys
import glob
import re
from optparse import OptionParser
from flair.Data import *


version = "1.2"


parser = OptionParser()
parser.add_option("-s", "--suffix", dest="suffix",
                  help="add this suffix to generated data. Default is \"_ascii.dat\".",
                  metavar="string")
parser.add_option("-c", "--clean", action="store_true", dest="clean",
                  default=False, 
                  help="Remove all files containing the suffix. Use carefully.")
parser.add_option("-d", "--detector", dest="detector",
                  help="If more detectors in binary file present, convert only this one.",
                  metavar="int")
parser.add_option("-F", "--FLUKA", action="store_true", dest="FLUKA",
                  default=False, 
                  help="Fluka style fortran output, 10 data points per line.")
(options, args) = parser.parse_args()

if options.detector is not None:
    detector_select = int(options.detector)
else:
    detector_select = -1

if len(args) < 1:
    print("")
    print("usrbin2ascii.py Version", version)
    print("")

    print("Usage: usrbin2ascii.py [options] binaryforfile")
    print("See: -h for all options.")
    print("File lists are supported when using quotations:")
    print("usrbin2ascii.py \"foo*\"")
    raise IOError("Error: no input file(s) stated.")

filename = args[0]
filename_list = glob.glob(filename)

if options.suffix is None:
    suffix = "_ascii.dat"
else:
    suffix = options.suffix


# cleanup all files with suffix
if options.clean:
        filename_list_suffix = glob.glob(filename+suffix)
        for filename_temp in filename_list_suffix:
                print("Removing filename", filename_temp)
                os.remove( filename_temp )


# reread filename list
# filename_list = glob.glob(filename)
# print filename_list

if len(filename_list) < 1:
    raise IOError("Error: %s does not exist." % filename)       

for filename in filename_list:
    print("opening binary file:", filename)
    mymatch = re.search(suffix,filename)
    if mymatch is not None:
        print()
        print("Hmmm. It seems you have not deleted previous ascii output.")
        print("I will exit now. You may want to wish to delete all files ending with",suffix)
        print("You can do this easily by rerunning the program and using the -c option.")
        sys.exit(0)


    print("="*80)
    usr = Usrbin(filename)
    usr.say()  # file,title,time,weight,ncase,nbatch
    for i in range(len(usr.detector)):
        print("-"*20,"Detector number %i" %i,"-"*20)
        usr.say(i) # details for each detector
    data = usr.readData(0)


    print("len(data):", len(data))
    fdata = unpackArray(data)
    print("len(fdata):", len(fdata))

    if len([x for x in fdata if x>0.0]) > 0:                    
        fmin = min([x for x in fdata if x>0.0])
        print("Min=",fmin)
    else:
        print("How sad. Your data contains only zeros.")
        print("Converting anyway.")
    if len(fdata) > 0:
        fmax = max(fdata)
        print("Max=",fmax)
    print("="*80)

    # TODO: handle multiple detectors.
    outfile = filename+suffix
    print("Writing output to", outfile)
    f = open(outfile,'w')
    det_number = 0
    
    pos = 0 # counter for data
    # TODO, check with multiple detectors. Are the data sequentially ordered?
    # At least we assume this here.

    if (detector_select + 1) > len(usr.detector):
        raise IOError("Selected detector number %i is not available." %detector_select)
    
    for det in usr.detector:
        #    det = usr.detector[0]
        # this header is only needed, if more detectors are requested
        if (len(usr.detector) > 1) and (detector_select == -1):
            f.write("# Detector number %i\n" % det_number )

        print("nx,ny,nz:", int(det.nx), int(det.ny), int(det.nz))

        ifort = 0 # counter for fortran format
        for iz in range(int(det.nz)):
            for iy in range(int(det.ny)):
                for ix in range(int(det.nx)):
                    fx = (ix/float(det.nx) * (det.xhigh-det.xlow)) + det.xlow
                    fy = (iy/float(det.ny) * (det.yhigh-det.ylow)) + det.ylow
                    fz = (iz/float(det.nz) * (det.zhigh-det.zlow)) + det.zlow
                    # middle of bins
                    fx += det.dx /2.0
                    fy += det.dy /2.0
                    fz += det.dz /2.0
                    if options.FLUKA:
                        output_string = "%.3e " % (fdata[pos])
                        ifort +=1
                        if ifort >= 10:
                            ifort = 0
                            output_string += "\n" # CR after 10 data points
                    else:
                        output_string = "%e %e %e %e\n" % (fx,fy,fz,fdata[pos])
                    pos += 1
                    if (det_number == detector_select) or (detector_select == -1):
                        f.write(output_string)
        det_number += 1
        f.write("\n")
    f.close()
