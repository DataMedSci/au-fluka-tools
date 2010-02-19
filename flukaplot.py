# -*- coding: utf-8 -*-
#!/usr/bin/env python
#########################################################################
# Copyright (C) 2010 David Hansen - davidchansen@hotmail.com
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
# 
# Plotting functions for FLUKA usrbins. 
# Note: For use with ipython, you must call ipython with -pylab
#
#
# 
#
# 
#
#
# TODO:
# - import the Data.py from flair in a better way
# - Document code better
import sys
import subprocess
import os
import math
import struct
from pylab import *

from scipy.io.numpyio import fread, fwrite
from numpy import *
from matplotlib.colors import LogNorm

from enthought.mayavi import mlab

#Warning, potentially very ugly way of importing the flair libs. 
flair = subprocess.Popen(["whereis","flair"],stdout=subprocess.PIPE).stdout.read().split(" ")[2]
flair = flair[0:len(flair)-1]
sys.path.append(flair)
sys.path.append(flair+"/lib")
    
from Data import *    

      
class AUUsrbin(Usrbin):
	
    
    	def arrayRead(self,det):
		
		bin = self.detector[det]
		data = self.readData(det)
		data  = asarray(unpackArray(data)).reshape((bin.nx,bin.ny,bin.nz), order='F')
		#Add data to the detector
		bin.data =data
		bin.start = [bin.xlow,bin.ylow,bin.zlow]
		bin.stop = [bin.xhigh,bin.yhigh,bin.zhigh]
		bin.delta = [bin.dx,bin.dy,bin.dz]
		bin.num = [bin.nx,bin.ny,bin.nz]
		
		#Compensate for normalization when not region binning
		if (bin.type in [1,7,11,17]): #Cylinder bins
			r = linspace(bin.xlow,bin.xhigh, num=bin.nx, endpoint=False)
			normalization = (pi*((2*bin.dx)*r+bin.dx**2)*bin.dz*bin.dy/(2*pi)).reshape(bin.nx,1,1)
			data = data*normalization 
		elif (bin.type in [0,3,4,5,10,13,14,15,16]): #Cartesian bins 
			normalization = bin.dx*bin.dy*bin.dz
			data = data*normalization 
		
		return data
		
			
	def plot(self,det,axis,start=None,stop=None):
		
		bin = self.detector[det]
		try:
			bin.data
		except AttributeError:
			self.arrayRead(det)
		sumaxis = [2,1,0]
		sumaxis.remove(axis)
		low,high = self.__ranges(det,start,stop)
		#Get requested data
		subdata = (slice(low[0],high[0]),slice(low[1],high[1]),slice(low[2],high[2]))
		plotdata = bin.data[subdata]
			
		for x in sumaxis:
			plotdata = sum(plotdata,x)
			
		#Renormalize data (renormalization theory not required)
		print str(low) + " " +str(high)
		if (bin.type in [1,7,11,17]): #Cylinder bins
			plotdata = plotdata/((high[0]**2*high[2]*high[1]/(2*pi))-(low[0]**2*low[2]*low[1]/(2*pi)))
		elif (bin.type in [0,3,4,5,10,13,14,15,16]): #Cartesian bins 
			plotdata = plotdata/((high[0]-low[0])*(high[1]-low[1])*(high[1]-low[1]))
		
		xdata = linspace(bin.start[axis]+bin.delta[axis]/2,bin.stop[axis]+bin.delta[axis]/2, num=bin.num[axis], endpoint=False)[slice(low[axis],high[axis])]
		plot(xdata,plotdata,label=bin.name)
		
	def plot2D(self,det,axisX,axisY,start=None,stop=None,logz=False):
		bin = self.detector[det]
		try:
			bin.data
		except AttributeError:
			self.arrayRead(det)
		sumaxis = [2,1,0]
		sumaxis.remove(axisX)
		sumaxis.remove(axisY)
		low,high = self.__ranges(det,start,stop)
		subdata = (slice(low[0],high[0]),slice(low[1],high[1]),slice(low[2],high[2]))
		plotdata = bin.data[subdata]
		plotdata = sum(plotdata,sumaxis[0])
		#Renormalize data
		if (bin.type in [1,7,11,17]): #Cylinder bins
			plotdata = plotdata/((high[0]**2*high[2]*high[1]/(2*pi))-(low[0]**2*low[2]*low[1]/(2*pi)))
		elif (bin.type in [0,3,4,5,10,13,14,15,16]): #Cartesian bins 
			plotdata = plotdata/((high[0]-low[0])*(high[1]-low[1])*(high[1]-low[1]))
			
		xdata = linspace(bin.start[axisX]+bin.delta[axisX]/2,bin.stop[axisX]+bin.delta[axisX]/2, num=bin.num[axisX], endpoint=False)[slice(low[axisX],high[axisX])]
		ydata = linspace(bin.start[axisY]+bin.delta[axisY]/2,bin.stop[axisY]+bin.delta[axisY]/2, num=bin.num[axisY], endpoint=False)[slice(low[axisY],high[axisY])]
		X,Y = meshgrid(ydata,xdata)
		if log:
			pcolormesh(X,Y,plotdata,norm=LogNorm(vmin=plotdata[plotdata.nonzero()].min(), vmax=plotdata.max()))
		else:
			pcolormesh(X,Y,plotdata)
	#Plots 2D data 	
	def plotSurf(self,det,axisX,axisY,start=None,stop=None):
		
		bin = self.detector[det]
		try:
			bin.data
		except AttributeError:
			self.arrayRead(det)
		sumaxis = [2,1,0]
		sumaxis.remove(axisX)
		sumaxis.remove(axisY)
		low,high = self.__ranges(det,start,stop)
		subdata = (slice(low[0],high[0]),slice(low[1],high[1]),slice(low[2],high[2]))
		plotdata = bin.data[subdata]
		plotdata = sum(plotdata,sumaxis[0])
		#Renormalize data
		if (bin.type in [1,7,11,17]): #Cylinder bins
			plotdata = plotdata/((high[0]**2*high[2]*high[1]/(2*pi))-(low[0]**2*low[2]*low[1]/(2*pi)))
		elif (bin.type in [0,3,4,5,10,13,14,15,16]): #Cartesian bins 
			plotdata = plotdata/((high[0]-low[0])*(high[1]-low[1])*(high[1]-low[1]))
			
		x = linspace(bin.start[axisX]+bin.delta[axisX]/2,bin.stop[axisX]+bin.delta[axisX]/2, num=bin.num[axisX], endpoint=False)[slice(low[axisX],high[axisX])]
		y = linspace(bin.start[axisY]+bin.delta[axisY]/2,bin.stop[axisY]+bin.delta[axisY]/2, num=bin.num[axisY], endpoint=False)[slice(low[axisY],high[axisY])]
		mlab.figure(1, fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
		z =plotdata
		# Visualize the points
		mlab.surf(x,y,z)
		mlab.show()
		
	
	def __ranges(self,det,start=None,stop=None):
		bin = self.detector[det]
		if start==None : start = bin.start
		if stop ==None : stop = bin.stop
		low = map(floor,[(x-xlow)/d for x,xlow,d in zip(start,bin.start,bin.delta)])
		low = [int(max(x,0)) for x in low]
		high = map(ceil,[(x-xlow)/d for x,xlow,d in zip(stop,bin.start,bin.delta)])
		high = [int(min(x,y)) for x,y in zip(high,bin.num)]
		return low,high