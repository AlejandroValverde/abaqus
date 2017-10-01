import os
import sys
import getopt

def recursiveFunction():
	cwd = os.getcwd() #Get working directory
	print('-> Exploring: '+cwd)
	for file in os.listdir(cwd):

		if os.path.isdir(file) and not 'git' in file:
			os.chdir(cwd + '/' + file)
			recursiveFunction()
			os.chdir(cwd)

		elif file.startswith('Thumbs.db'):
			print('--> Thumbs removed in: '+cwd)
			os.remove(file)

		elif file.endswith('.py') or file.endswith('.txt') or file.endswith('.dat'):# and not 'eps' in cwd:

			print('--> Converting: '+file+ ', in: '+cwd)
			os.system('dos2unix '+file)
			# if not os.path.isdir('eps'):
			# 	os.mkdir('.\\eps')
				
			# print('--> Moving file: '+file)
			# try:
			# 	os.rename(cwd+'\\'+file, cwd+'\\eps\\'+file)
			# except FileExistsError as e:
			# 	print('--> Overwriting file: '+file)
			# 	os.remove(cwd+'\\eps\\'+file)
			# 	os.rename(cwd+'\\'+file, cwd+'\\eps\\'+file)


recursiveFunction()

print('---> Execution finished')