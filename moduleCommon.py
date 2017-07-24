import os
import platform
import pdb #pdb.set_trace()
from shutil import copyfile

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