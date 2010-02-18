#!/usr/bin/env python
# $Id: Data.py,v 1.11 2007/05/10 16:15:13 bnv Exp $
#
# Copyright and User License
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
# Copyright Vasilis.Vlachoudis@cern.ch for the
# European Organization for Nuclear Research (CERN)
#
# All rights not expressly granted under this license are reserved.
#
# Installation, use, reproduction, display of the
# software ("flair"), in source and binary forms, are
# permitted free of charge on a non-exclusive basis for
# internal scientific, non-commercial and non-weapon-related
# use by non-profit organizations only.
#
# For commercial use of the software a license fee should be
# payed to the author. Please contact the main author
# Vasilis.Vlachoudis@cern.ch for further information.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following
# conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
#
# DISCLAIMER
# ~~~~~~~~~~
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, IMPLIED WARRANTIES OF MERCHANTABILITY, OF
# SATISFACTORY QUALITY, AND FITNESS FOR A PARTICULAR PURPOSE
# OR USE ARE DISCLAIMED. THE COPYRIGHT HOLDERS AND THE
# AUTHORS MAKE NO REPRESENTATION THAT THE SOFTWARE AND
# MODIFICATIONS THEREOF, WILL NOT INFRINGE ANY PATENT,
#!/usr/bin/env python
#########################################################################
# Copyright (C) 2009 Niels Bassler, bassler@phys.au.dk
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
# change log
# 18.2.2010, NB
# - removed flair code, which is now imported instead.
# 22.11.2008, NB
# - added fortran block length 128 support from upstream release
# - added -F option
# 23.9.2008, NB
# - added help text and option parser
# - output suffix can now be set optionally and cleanup option
# - added glob.glob
# - fix exception if datafiles only contains zeros
#
#
# TODO:
# - check if input files are binary before processing, improve file checking
# - allow usrbin2ascii foo*   instead of     usrbin2ascii "foo*"
# - less verbose
# - merge usrbin2ascii wtih usrtrack2ascii
# - convert to shell scripts which sets flair path for python instead of this solution ?



import math
import struct

#NB
import os,sys, glob,re
from optparse import OptionParser

def get_datapy_path():
    flairexe = "flair"
    datapy = "Data.py"
    path = os.environ["PATH"]
    paths = path.split(os.pathsep)
    extlist = ['']
    for p in paths:
        f = os.path.join(p, flairexe)
	d = os.path.join(p, datapy)
	if os.path.isfile(f):
            if os.path.isfile(d):
                return(p)
	    else:
                raise IOError("Could not find Data.py")
	
    raise IOError("Could not find flair installation in PATH.")


# ----------------------------------------------
if __name__ == "__main__":
	ddd = get_datapy_path()
	sys.path.append(ddd)
	sys.path.append(ddd+"/lib")
	print sys.path
	
	from Data import *

	version = "1.2"

	parser = OptionParser()
	parser.add_option("-s", "--suffix", dest="suffix",
                  help="add this suffix to generated data. Default is \"_ascii.dat\".", metavar="string")
	parser.add_option("-c", "--clean", action="store_true", dest="clean",
                  default=False, help="Remove all files containing the suffix. Use carefully.")
	parser.add_option("-F", "--FLUKA", action="store_true", dest="FLUKA",
                  default=False, help="Fluka style fortran output, 10 data points per line.")


	(options, args) = parser.parse_args()

	if len(args) < 1:
		print ""
		print "usrbin2ascii.py Version", version
		print ""
		print "Error: no input file stated."
		print "Usage: usrbin2ascii.py [options] binaryforfile"
		print "See: -h for all options."
		print "Filelists are supported when using quotations:"
		print "usrbin2ascii.py \"foo*\""
		sys.exit(0)
	
	filename = args[0]

	if options.suffix == None:
		suffix = "_ascii.dat"
	else:
		suffix = options.suffix
#	print suffix
#	print args[0]

	filename_list = glob.glob(filename)

	# cleanup all files with suffix
	if options.clean == True:
		filename_list = glob.glob(filename+suffix)
		for filename_temp in filename_list:
			print "Removing filename", filename_temp
			os.remove( filename_temp )




	# reread filename list
	filename_list = glob.glob(filename)
	print filename_list
	
	if len(filename_list) < 1:
		print("Error: %s does not exist." % filename)	
		sys.exit(0)
		

#	print filename, filename_list, len(filename_list)





	for filename in filename_list:
		print "opening binary file:", filename
		mymatch = re.search(suffix,filename)
		if mymatch != None:
			print
			print "Hmmm. It seems you have not deleted previous ascii output."
			print "I will exit now. You may want to wish to delete all files ending with", suffix
			print "You can do this easily by rerunning the program and using the -c option."
			sys.exit(0)
			
		


		print "="*80
		usr = Usrbin(filename)
		usr.say()
		for i in range(len(usr.detector)):
			print "-"*50
			usr.say(i)


		data = usr.readData(0)
	# no statistics there.
	#	stats = usr.readStats(0)

		print "len(data):", len(data)
		fdata = unpackArray(data)
		print "len(fdata):", len(fdata)

		if len([x for x in fdata if x>0.0]) > 0:			
			fmin = min([x for x in fdata if x>0.0])
			print "Min=",fmin
		else:
			print "How sad. Your data contains only zeros."
			print "Converting anyway."
		if len(fdata) > 0:
			fmax = max(fdata)
			print "Max=",fmax
	# 	print "Tot=",total[0],total[1]


		print "="*80

		det = usr.detector[0]	
		pos = 0
		
		outfile = filename+suffix

		print "nx,ny,nz:", int(det.nx), int(det.ny), int(det.nz)

		print "Writing output to", outfile
		f = open(outfile,'w')
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
					#debug
					#output_string = "%i %i %i %i %e %e %e %e\n" % (pos, ix,iy,iz,fx,fy,fz,fdata[pos])
					if options.FLUKA == True:
						output_string = "%.3e " % (fdata[pos])
						ifort +=1
						if ifort >= 10:
							ifort = 0
							output_string += "\n" # CR after 10 data points
					else:
						output_string = "%e %e %e %e\n" % (fx,fy,fz,fdata[pos])
					pos += 1
					f.write(output_string)
		f.close()


