import os
import platform
import pdb #pdb.set_trace()
from shutil import copyfile
import numpy as np

def ConvertNumber(string):

	string = string.replace('\n','')

	string = string.lstrip()

	string = string.rstrip()

	index = string.find('E')

	if index == -1: #if number does not contain "E"

		return float(string)

	else:

		num = float(string[:index])

		exp = float(string[(index+1):])

		return num * 10**(exp)

def getQandPvectors(design):

	totalLength = design.cutWingTip - design.cutWingRoot

	totalHeight = design.cutUp + abs(design.cutDown)

	#The vector Q iterates through the first set of columns of chiral nodes
	Q_i = np.arange(design.distanceCenterPoints, totalLength + design.distanceCenterPoints, design.distanceCenterPoints)
	Q_j = np.arange(0.0, totalHeight, 2 * design.heightTriangle)

	#The vector P iterates through the second set of columns of chiral nodes
	P_i = np.arange(design.distanceCenterPoints * 3/2, totalLength + design.distanceCenterPoints, design.distanceCenterPoints)
	P_j = np.arange(design.heightTriangle, totalHeight, 2 * design.heightTriangle)

	return Q_i, Q_j, P_i, P_j

def globalChangeDir(cwd, address):

	if address == '.':
		os.chdir(cwd)

	else:

		if platform.system() == 'Linux':
		    newAddress = address.replace('-', '/')      
		elif platform.system() == 'Windows':
		    newAddress = address.replace('-', '\\')
		else:
		    newAddress = address.replace('-', '/')
		    print('OS not recognized, assumed unix based')

		os.chdir(cwd + newAddress)

def globalCreateDir(cwd, address):

	if platform.system() == 'Linux':
	    newAddress = address.replace('-', '/')      
	elif platform.system() == 'Windows':
	    newAddress = address.replace('-', '\\')
	else:
	    newAddress = address.replace('-', '/')
	    print('OS not recognized, assumed unix based')
	if not os.path.isdir('.'+newAddress):
		os.mkdir('.'+newAddress)

def globalCopyFile(address1, address2, fileName1, fileName2):

	if platform.system() == 'Linux':
		link = '/'
		newAddress1 = address1.replace('-', '/') 
		newAddress2 = address2.replace('-', '/') 
	elif platform.system() == 'Windows':
		link = '\\'
		newAddress1 = address1.replace('-', '\\') 
		newAddress2 = address2.replace('-', '\\') 
	else:
		link = '/'
		newAddress1 = address1.replace('-', '/') 
		newAddress2 = address2.replace('-', '/') 
		print('OS not recognized, assumed unix based')

	copyfile(newAddress1+link+fileName1, newAddress2+link+fileName2)

def isUnix():

	if platform.system() == 'Linux':
	    return True    
	elif platform.system() == 'Windows':
	    return False
	else:
	    return True
	    print('OS not recognized, assumed unix based')