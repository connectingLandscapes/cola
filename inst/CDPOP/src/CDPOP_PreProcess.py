# -------------------------------------------------------------------------------------------------
# CDPOP_PreProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file pre processing.
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError("Numpy required.")

# Scipy function KDTree
try:
	from scipy.spatial import KDTree
	scipyAvail = True
except ImportError:
	raise ImportError("Scipy required.")
	scipyAvail = False
	
# CDPOP functions
try:
	from CDPOP_Mate import *
	from CDPOP_Offspring import *
	from CDPOP_Disperse import *
except ImportError:
	raise ImportError("CDPOP Modules missing.")
	
# General python imports
import os,sys

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ----------------------------------------------------------------------------------
def loadFile(filename, header_lines=0, delimiter=None, cdpop_inputvars=False): ###
	'''
	Used to load file hearders according to current UNICOR standards
	as of 5/19/2011
	'''
	try:
		inputfile = open(filename)
	except (IOError,OSError) as e:
		print(("Load file: %s the file (%s) is not available!"%(e,filename)))
		sys.exit(-1)
	header_dict = {}
	data_list = []
	index_list = [] ###
	
	if delimiter != None:
		lines = [ln.rstrip().split(delimiter) for ln in inputfile]
	else:
		lines = [ln.rstrip().split() for ln in inputfile]
	# Close file
	inputfile.close()
		
	for i,line in enumerate(lines):
		if i < header_lines:
			if len(line) <= 1:
				print("Only one header value in line, skipping line...")
				continue
			#else:	
			elif len(line) == 2: ###
				header_dict[line[0]] = line[1]
			### working with a file where the first line is all header keys
			### and the following lines are their data
			elif cdpop_inputvars: ###
				for j in range(len(line)): ###
					header_dict[line[j]] = [] ###
					index_list.append(line[j]) ###
		else:
			#This is a fix to remove empty entries from from a line if they
			#exist, this is to fix a problem with reading in cdmatrices
			for i in range(line.count('')):
				line.remove('')
			data_list.append(line)
			if cdpop_inputvars: ###
				#tempTuple = ()
				for j in range(len(line)): ###
					# remove the following lines should the tuple representing bar-delimited values break anything -TJJ
					if line[j].find('|') != -1:
						tempList = line[j].split('|')
						line[j] = tuple(tempList)
					#---
					
					header_dict[index_list[j]].append(line[j]) ###
	
	if not cdpop_inputvars:
		return header_dict, data_list
	else:
		n_jobs = len(lines) - header_lines
		return header_dict, index_list, n_jobs

# ---------------------------------------------------------------------------------------------------	 
def GetMaxCDValue(threshold,cdmatrix):	
	'''
	GetMaxCDValue()
	This function calculates the maximum movement thresholds.
	'''	
	
	# movement threshold if max specified
	if str(threshold).endswith('max'):
		# If max
		if len(threshold.strip('max')) == 0:
			threshold = np.amax(cdmatrix)
		else:
			threshold = (float(threshold.strip('max'))/100.)\
			*np.amax(cdmatrix)
	else:
		threshold = float(threshold)
	
	return threshold
# ---------------------------------------------------------------------------------------------------	 
def ReadFitnessSurface(fitsurface):	
	'''
	ReadFitnessSurface()
	This function reads in the ascii fitness surface.
	'''	
	# Open file for reading
	inputfile = open(fitsurface,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	values = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		i = i.strip('\n').strip('\r').strip(' ')
		thisline = i.split(' ')
		values.append(thisline)
	
	# Grab some information from fitvalues: number of columns
	lenfirstline = len(values[0])
	ncols = int(values[0][lenfirstline-1])
		
	# Grab some information from values: number of rows
	lensecondline = len(values[1])
	nrows = int(values[1][lensecondline-1])

	# Grab some information from values: x coord value in lower left
	lenthirdline = len(values[2])
	xllcorner = float(values[2][lenthirdline-1])

	# Grab some information from values: y coord value in lower left
	lenfourthline = len(values[3])
	yllcorner = float(values[3][lenfourthline-1])

	# Grab some information from values: cell size
	lenfifthline =len(values[4])
	cellsize = float(values[4][lenfifthline-1])
	
	# Grab some information from values: Nodataval
	lensixthline =len(values[5])
	Nodataval = float(values[5][lensixthline-1])
	
	#xll,yll is the bottom left corner location, want key spot locations to be the center of each cell
	xll = xllcorner + (cellsize/2.)
	yll = yllcorner + (cellsize/2.)
	
	# Error statement
	if len(values)-6 != nrows or len(values[6]) != ncols:
		print('Spatial selection surface file is not in correct format. Return error number of dimensions.')
		sys.exit(-1)
	# And turn rasterfile into a list with float values without header information
	fitnessvalues = []
	for i in range(len(values)-6):
		fitnessvalues.append([])
		for j in range(len(values[6])):
			fitnessvalues[i].append(float(values[6+i][j]))
	
	# Return tuple
	tupReadFit = fitnessvalues,ncols,nrows,xll,yll,cellsize,Nodataval
	return tupReadFit
	# End::ReadFitnessSurface()
	
# ---------------------------------------------------------------------------------------------------	
def CreateAlleleList(loci,alleles,xgenes):
	'''
	CreateAlleleList()
	This function creates a list for the allele storage.
	'''	
	
	# Store all information in a list [loci][allele#,probability]
	allelst = []
	for i in range(loci):
		allelst.append([])
		for k in range(alleles[i]):
			#allspot = alleles[i]*i+1+k
			allspot = sum(alleles[0:i]) + k + 1
			allelst[i].append([int(k),float(xgenes[allspot][1])])
	
	# Return variables
	return allelst
	
	# End::CreateAlleleList()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeGenes(intgenesans,allefreqfilename,loci,alleles,datadir):	
	
	allelst = []
	
	# If genetic structure intialized by a file...
	if intgenesans == 'file' or intgenesans == 'file_introduce' or intgenesans == 'file_var' or intgenesans == 'file_introduce_var':
		
		# Loop through allele frequency files given
		for ifile in range(len(allefreqfilename)):
			
			fileans = allefreqfilename[ifile]
			
			# If genetic structure intialized by a file...
			if fileans != 'random':
				
				# Check statements
				if os.path.exists(datadir+fileans+'.csv'):
					# Open file for reading
					inputfile = open(datadir+fileans+'.csv','r')
				else:
					print(("CDPOP InitializeGenes() error: open failed, could not open %s"%(fileans)))
					sys.exit(-1)
						
				# Read lines from the file
				lines = inputfile.readlines()
				
				#Close the file
				inputfile.close()
				
				# Create an empty matrix to append to
				xgenes = []
			
				# Split up each line in file and append to empty matrix, x
				for i in lines:
					thisline = i.strip('\n').strip('\r').strip(' ').split(',')
					xgenes.append(thisline)
				
				# Error check here
				if (len(xgenes)-1) != sum(alleles):
					print('Allele frequency file is not the specified number of loci and alleles as in in input file.')
					sys.exit(-1)
				
				# Delete lines from earlier
				del(lines)
				
				# Call CreateAlleleList()
				allelst.append(CreateAlleleList(loci,alleles,xgenes))
			
			# If genetic structure is to be initialize by random
			elif fileans == 'random':
				
				# Create even distribution
				xgenes = []
				xgenes.append(['Allele List','Frequency'])
				for iloci in range(loci):
					for iall in range(alleles[iloci]):
						xgenes.append(['L'+str(iloci)+'A'+str(iall),str(1.0/alleles[iloci])])
				
				# Call CreateAlleleList()
				allelst.append(CreateAlleleList(loci,alleles,xgenes))	
	
	# If just random for entire population
	elif intgenesans == 'random' or intgenesans == 'random_var':
		
		# Create even distribution
		xgenes = []
		xgenes.append(['Allele List','Frequency'])
		for iloci in range(loci):
			for iall in range(alleles[iloci]):
				xgenes.append(['L'+str(iloci)+'A'+str(iall),str(1.0/alleles[iloci])])
		
		# Call CreateAlleleList()
		allelst.append(CreateAlleleList(loci,alleles,xgenes))
	
		# Delete x variable
		del(xgenes)
	
	# Error statement
	else:
		allelst = []
			
	# Return variables
	return allelst
		
	# End::InitializeGenes()
	
# ---------------------------------------------------------------------------------------------------	 
def InitializeAge(agefilename,nogrids,datadir):
	'''
	InitializeAge()
	This function initializes the age of each population
	with an age distribution list.
	'''
	
	# Loop through age files given
	agelst = []
	Magemort = []
	Fagemort = []
	eggs_mean = []
	eggs_sigma = []	
	Mmature = []
	Fmature = []
	Mnewmortperc = []
	Fnewmortperc = []
	for ifile in range(len(agefilename)):
		agelst.append([])
		Magemort.append([])
		Fagemort.append([])
		eggs_mean.append([])
		eggs_sigma.append([])	
		Mmature.append([])
		Fmature.append([])
		Mnewmortperc.append([])
		Fnewmortperc.append([])
		
		fileans = agefilename[ifile]			
				
		# Check statements
		if os.path.exists(datadir+fileans):
			# Open file for reading
			inputfile = open(datadir+fileans,'r')
		else:
			print(("CDPOP InitializeAge() error: open failed, could not open %s"%(datadir+fileans)))
			sys.exit(-1)
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Error statement on agevars file length
		if len(lines) < 3:
			print('Agevars.csv file must have at least 2 age classes: Age 0 (eggs) and Age 1 (Adults).')
			sys.exit(-1)
		
		# Create an empty matrix to append to
		xage = []
		
		# Split up each line in file and append to empty matrix, x
		for i in lines:
			thisline = i.strip('\n').strip('\r').strip(' ').split('\r')
			for j in range(len(thisline)):
				xage.append(thisline[j].split(','))
			
		# Delete lines from earlier
		del(lines)
		
		# Store all information in a list [age,probability]
		ageclass = []
		ageno = []
		for i in range(len(xage)-2): # i+2, only get Age 1+
			ageclass.append(int(xage[i+2][0]))
			ageno.append(float(xage[i+2][1]))
			if len(xage[i+2][2].split('|')) == 1: 
				Magemort[ifile].append(float(xage[i+2][2])/100.)
			else:
				Magemort[ifile].append(xage[i+2][2].split('|'))
			if len(xage[i+2][3].split('|')) == 1: 
				Fagemort[ifile].append(float(xage[i+2][3])/100.)
			else:
				Fagemort[ifile].append(xage[i+2][3].split('|'))
			if len(xage[i+2][4].split('|')) == 1: 
				eggs_mean[ifile].append(float(xage[i+2][4]))
			else:
				eggs_mean[ifile].append(xage[i+2][4].split('|'))
			if len(xage[i+2][5].split('|')) == 1: 
				eggs_sigma[ifile].append(float(xage[i+2][5]))
			else:
				eggs_sigma[ifile].append(xage[i+2][5].split('|'))
			if len(xage[i+2][6].split('|')) == 1:
				Mmature[ifile].append(float(xage[i+2][6]))
			else:
				Mmature[ifile].append(xage[i+2][6].split('|'))
			if len(xage[i+2][7].split('|')) == 1:
				Fmature[ifile].append(float(xage[i+2][7]))
			else:
				Fmature[ifile].append(xage[i+2][7].split('|'))
		
		# Get age distribution list
		for i in range(len(ageno)):
			agelst[ifile].append([ageclass[i],ageno[i]/sum(ageno)])
		
		# Get Age 0 mortality for males and females
		if len(xage[1][2].split('|')) == 1: 
			Mnewmortperc[ifile] = float(xage[1][2])/100.
		else:
			Mnewmortperc[ifile] = xage[1][2].split('|')
		if len(xage[1][3].split('|')) == 1: 
			Fnewmortperc[ifile] = float(xage[1][3])/100.
		else:
			Fnewmortperc[ifile] = xage[1][3].split('|')
			
		# Error checks here: if number of classes does not equal mortality age classes
		if len(agelst[ifile]) != len(Magemort[ifile]):
			print('Agedistribution data not fully entered correctly.')
			sys.exit(-1)
		# Error check that age 0s are initialized
		if xage[1][1] != '0':
			print('Warning: Age 0 was initialized at start of model, enter 0 for distribution in Agevars.csv file.')
			sys.exit(-1)
		# Deletes
		del(xage)
	
	# Return variables
	return agelst,Magemort,Fagemort,eggs_mean,eggs_sigma,Mnewmortperc,Fnewmortperc,Mmature,Fmature
	
	# End::InitializeAge()

# ---------------------------------------------------------------------------------------------------	 
def ReadXY(xyfilename):
	'''
	ReadMateXYCDMatrix()
	This function reads in the xy values for the cost distance matrix.
	'''	
	
	
	# Open file for reading
	inputfile = open(xyfilename+'.csv','r')
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	#Close the file
	inputfile.close()
	
	# Create an empty matrix to append to
	xy = []
	
	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.rstrip('\n').rstrip('\r').split(',')
		xy.append(thisline)
		
	# Delete lines from earlier
	del(lines)
	
	# Return variables
	return xy
	
	# End::ReadXY()

# ---------------------------------------------------------------------------------------------------	 
def ReadCDMatrix(cdmatrixfilename,function,threshold,A,B,C,subpop):
	'''
	ReadMateCDMatrix()
	This function reads in the mating cost distance matrix.
	'''	
	
	# Check statements
	if os.path.exists(cdmatrixfilename+'.csv'):
		# Open file for reading
		inputfile = open(cdmatrixfilename+'.csv','r')
	else:
		print(("CDPOP ReadCDMatrix() error: open failed, could not open %s"%(cdmatrixfilename+'.csv')))
		sys.exit(-1)
	
	# Read lines from the file
	lines = inputfile.readlines()
	
	# Close the file
	inputfile.close()
	
	# Create an empty matrix to append to 
	bigCD = []
	
	# Split up each line in file and append to empty matrix, x
	for spot in lines:
		thisline = spot.strip('\n').split(',')
		bigCD.append(thisline[0:len(lines)])
	bigCD = np.asarray(bigCD,dtype='float')
	
	# Delete lines from earlier
	del(lines)
		
	# Store number of files
	nofiles = len(bigCD)
	
	# Calculate max and min of bigCD matrix
	minbigCD = np.amin(bigCD)
	maxbigCD = np.amax(bigCD)
	
	# Error checks on these values
	if minbigCD < 0:
		print('Cost matrix values should not be negative.')
		sys.exit(-1)
	if maxbigCD != 0:
		if maxbigCD < 1 and (function == '2' or function == '5'):
			print('Careful use of cost distance values less than 1 when using options 2 or 5.')
			sys.exit(-1)
			
	# Get maximum cdvalue to use for movethreshold if specified
	threshold = GetMaxCDValue(threshold,bigCD)
	
	# Create a matrix of to be filled 
	cdmatrix = []
	
	# If not the prob matrix directly
	if function != '9':
		# Fill up matrix with float value of array x
		for j in range(nofiles):
			cdmatrix.append([])
			for k in range(nofiles):
				
				# For the linear function
				if function == '1':
					scale_min = 0.
					scale_max = threshold
					# Set = to 0 if cdvalue is greater than or equal to movethreshold
					if float(bigCD[j][k]) >= threshold:
						cdmatrix[j].append(0.0)
					# If threshold is 0 (philopatry) set to 1 - can't dived by 0
					elif float(bigCD[j][k]) < threshold and threshold == 0.0:
						cdmatrix[j].append(1.0)
					# Else calculated function value and if not philopatry
					elif float(bigCD[j][k]) < threshold and threshold != 0.0:
						cdmatrix[j].append(-(1./threshold)*float(bigCD[j][k]) + 1)
					else:
						print('Something off in linear function values.')
						sys.exit(-1)
							
				# For the inverse square function
				elif function == '2':
					
					# This function gets rescale: calculate here
					if threshold == 0:
						scale_min = 0.
					else:
						scale_min = 1./(pow(threshold,2))
					scale_max = 1.
					
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) >= threshold:
						cdmatrix[j].append(0.0)
					# If threshold is 0 (philopatry) set to 1 - can't dived by 0
					elif float(bigCD[j][k]) < threshold and threshold == 0.0:
						cdmatrix[j].append(1.0)
					# If cd mat is 0. 
					elif float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0 or (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						invsq_val = 1./(pow(float(bigCD[j][k]),2))
						#invsq_val = (invsq_val - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(invsq_val)# Else something else.
					else:
						print('Something off in inv squ function values.')
						sys.exit(-1)
						#probarray.append([fgspot,1./(pow(float(freeoffcd[fgspot]),2))])
								
								
				# Nearest neighbor function here
				elif function == '3':
					print('Nearest neighbor function is not currently implemented.')
					print('You can use Linear function with neighbor threshold for approximation. Email Erin.')
					sys.exit(-1)
					
				# Random function here
				elif function == '4':
					scale_min = 0.
					scale_max = threshold
					if threshold != 0: # If matrix is not all 0s
						# Set = to 0 if cdvalue is greater than movethreshold
						if float(bigCD[j][k]) >= threshold:
							cdmatrix[j].append(0.0)
						# Else calculated function value
						else:
							cdmatrix[j].append(1.0)
					else:
						cdmatrix[j].append(1.0)
				
				# For the negative binomial function
				elif function == '5':
				
					# This function gets rescale: calculate here
					scale_min = A*pow(10,-B*float(threshold))
					scale_max = A*pow(10,-B*float(minbigCD))
				
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) >= threshold:
						cdmatrix[j].append(0.0)
					# Rescaled value divide by zero check cases
					elif float(bigCD[j][k]) < threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) < threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						negexp = A*pow(10,-B*float(bigCD[j][k]))
						negexp = (negexp - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(negexp)
					# Else something else.
					else:
						print('Something off in neg exp function values.')
						sys.exit(-1)
						
				# For in a subpopulation only
				elif function == '6':
					
					# Skip if only one subpopulation
					if len(count_unique(np.asarray(subpop))[0]) != 1:
						# Check if within the same subpopulation
						if subpop[j] == subpop[k]:
							# This function gets rescale: calculate here
							if threshold == 0:
								scale_min = 0.
							else:
								scale_min = 1./(pow(threshold,2))
							scale_max = 1.
							
							# Set = to 0 if cdvalue is greater than movethreshold
							if float(bigCD[j][k]) >= threshold:
								cdmatrix[j].append(0.0)
							# If threshold is 0 (philopatry) set to 1 - can't dived by 0
							elif float(bigCD[j][k]) < threshold and threshold == 0.0:
								cdmatrix[j].append(1.0)
							# If cd mat is 0. 
							elif float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0 or (minbigCD == maxbigCD or int(maxbigCD) == 0):
								cdmatrix[j].append(1.0)
							# Else calculated function value
							elif float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
								invsq_val = 1./(pow(float(bigCD[j][k]),2))
								invsq_val = (invsq_val - scale_min) / (scale_max - scale_min)
								cdmatrix[j].append(invsq_val)# Else something else.
							else:
								print('Something off in inv squ function values.')
								sys.exit(-1)					
						else:
							cdmatrix[j].append(0.0)
					else:
						cdmatrix[j].append(1.0)
						scale_max = 1.
						scale_min = 0.
				
				# For Gaussian function 
				elif function == '7':
					if C == 0.0:
						print('Parameter C for Gaussian function is 0.')
						sys.exit(-1)
					# This function gets rescale: calculate here
					scale_min = A*np.exp(-((float(threshold)-B)**2)/(2*C**2))
					scale_max = A*np.exp(-((float(minbigCD)-B)**2)/(2*C**2))
				
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) >= threshold:
						cdmatrix[j].append(0.0)
					# Rescaled value divide by zero check cases
					elif float(bigCD[j][k]) < threshold and threshold == 0.0 and (minbigCD == maxbigCD or int(maxbigCD) == 0):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif float(bigCD[j][k]) < threshold and threshold != 0.0 and (minbigCD != maxbigCD and int(maxbigCD) != 0):
						gauss_val = A*np.exp(-((float(bigCD[j][k])-B)**2)/(2*C**2))
						gauss_val = (gauss_val - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(gauss_val)
					# Else something else.
					else:
						print('Something off in gauss function values.')
						sys.exit(-1)
						
				# For cost distance matrix only function 
				elif function == '8':
					
					scale_min = minbigCD
					scale_max = threshold
					
					# Set = to 0 if cdvalue is greater than movethreshold
					if float(bigCD[j][k]) >= threshold:
						cdmatrix[j].append(0.0) 
					# Rescaled value divide by zero check cases - philopatry
					elif (float(bigCD[j][k]) < threshold and threshold == 0.0) and (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
						cdmatrix[j].append(1.0)
					# If cd mat is 0. 
					elif (float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) == 0.0) or (minbigCD == maxbigCD or int(maxbigCD) == 0 or threshold == minbigCD):
						cdmatrix[j].append(1.0)
					# Else calculated function value
					elif (float(bigCD[j][k]) < threshold and threshold != 0.0 and float(bigCD[j][k]) != 0.0) and (minbigCD != maxbigCD and int(maxbigCD) != 0 and threshold != minbigCD):
						cd_val = (float(bigCD[j][k]) - scale_min) / (scale_max - scale_min)
						cdmatrix[j].append(1. - cd_val)
					# Else something else.
					else:
						print('Something off in 8 function values.')
						sys.exit(-1)
				
				# error
				else:
					print('This movement function option does not exist.')
					sys.exit(-1)
	else: # For option 9
		cdmatrix = bigCD
		scale_min = minbigCD
		scale_max = maxbigCD
	
	# Transpose all matrices here (this was because gdistance asymetrical matrices read by column
	cdmatrix = np.transpose(np.asarray(cdmatrix))
	
	# Delete variables
	del(bigCD)
	
	# Return variables
	tupReadMat = cdmatrix,threshold,scale_min,scale_max
	return tupReadMat
	
	# End::ReadMateCDMatrix

# ---------------------------------------------------------------------------------------------------	 
def DoCDClimate(datadir,icdtime,cdclimgentime,matecdmatfile,dispcdmatfile,matemoveno,Fdispmoveno,Mdispmoveno,matemovethresh,Fdispmovethresh,Mdispmovethresh,matemoveparA,matemoveparB,matemoveparC,FdispmoveparA,FdispmoveparB,FdispmoveparC,MdispmoveparA,MdispmoveparB,MdispmoveparC,subpop,Magemort,Fagemort,offno,lmbda,sigma,K_envvals,Mnewmortperc,Fnewmortperc,fitvals,twinning,Mmature,Fmature,betaFile_selection,xvars_betas_pass,epimod_pass,epireset_pass,betaFile_epigene,cdevolveans,epigeneans,gridmort_pass):
	'''
	DoCDClimate()
	Reads in cost distance matrices and converts to probabilities.
	'''
	
	# -------------------------------
	# Extract cdclimate values here
	# -------------------------------
	# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
	if not isinstance(cdclimgentime, (list,tuple)):
		# Error checks
		if cdclimgentime != ['0']:
			print('If not using CDClimate option, set begin time loop with cdclimgentime at 0.')
			sys.exit(-1)
	
	# --------------------------------------
	# For values linked to subpops (AgeVars)
	# --------------------------------------
	
	# Get age mortality here
	tempMmort = []
	tempFmort = []
	templmbda = []
	tempsigma = []
	tempMnewmort = []
	tempFnewmort = []
	tempMmat = []
	tempFmat = []
	# Loop through subpops
	for jpop in range(len(Magemort)):
		tempMmort.append([])
		tempFmort.append([])
		templmbda.append([])
		tempsigma.append([])
		tempMnewmort.append([])
		tempFnewmort.append([])
		tempMmat.append([])
		tempFmat.append([])
		# Loop through the ages
		for i in range(len(Magemort[jpop])):
			
			if not isinstance(Magemort[jpop][i], float):
				tempMmort[jpop].append(float(Magemort[jpop][i][icdtime])/100.)
			else:
				tempMmort[jpop].append(Magemort[jpop][i])
			if not isinstance(Fagemort[jpop][i], float):
				tempFmort[jpop].append(float(Fagemort[jpop][i][icdtime])/100.)
			else:
				tempFmort[jpop].append(Fagemort[jpop][i])
			if not isinstance(lmbda[jpop][i], float):
				templmbda[jpop].append(float(lmbda[jpop][i][icdtime]))
			else:
				templmbda[jpop].append(lmbda[jpop][i])
			if not isinstance(sigma[jpop][i], float):
				tempsigma[jpop].append(float(sigma[jpop][i][icdtime]))
			else:
				tempsigma[jpop].append(sigma[jpop][i])
			if not isinstance(Mmature[jpop][i], float):
				tempMmat[jpop].append(float(Mmature[jpop][i][icdtime]))
			else:
				tempMmat[jpop].append(Mmature[jpop][i])
			if not isinstance(Fmature[jpop][i], float):
				tempFmat[jpop].append(float(Fmature[jpop][i][icdtime]))
			else:
				tempFmat[jpop].append(Fmature[jpop][i])
	
		if not isinstance(Magemort[jpop][i], float):
			tempMnewmort[jpop] = float(Mnewmortperc[jpop][icdtime])/100.
		else:
			tempMnewmort[jpop] = Mnewmortperc[jpop]
		if not isinstance(Fagemort[jpop][i], float):
			tempFnewmort[jpop] = float(Fnewmortperc[jpop][icdtime])/100.
		else:
			tempFnewmort[jpop] = Fnewmortperc[jpop]
	
	# ----------------------
	# Patch based parameters
	# ----------------------
	tempfitvals = []
	tempxvars_betas = []
	tempepimod = []
	tempepireset = []
	tempgridmort = []
	for isub in range(len(subpop)):
		if len(fitvals) > 0:
			tempfitvals.append([])
			for i in range(len(fitvals[isub])):
				if len(fitvals[isub][i].split('|')) > 1:
					tempfitvals[isub].append(fitvals[isub][i].split('|')[icdtime])
				else:
					tempfitvals[isub].append(fitvals[isub][i])
				# Check if ; separated for 1_HetMort
				if len(tempfitvals[isub][i].split(';')) > 1:
					# Quick error check
					if len(tempfitvals[isub][i].split(';')) != 2:
						print('Two values used for 1_HeMort option.')
						sys.exit(-1)
					tempfitvals[isub][i] = tempfitvals[isub][i].split(';')
		if len(xvars_betas_pass) > 0:		
			tempxvars_betas.append([])
			for i in range(len(xvars_betas_pass[isub])):
				if len(xvars_betas_pass[isub][i].split('|')) > 1:
					tempxvars_betas[isub].append(xvars_betas_pass[isub][i].split('|')[icdtime])
				else:
					tempxvars_betas[isub].append(xvars_betas_pass[isub][i])
		if len(epimod_pass) > 0:		
			tempepimod.append([])
			for i in range(len(epimod_pass[isub])):
				if len(epimod_pass[isub][i].split('|')) > 1:
					tempepimod[isub].append(epimod_pass[isub][i].split('|')[icdtime])
				else:
					tempepimod[isub].append(epimod_pass[isub][i])	
		if len(epireset_pass) > 0:		
			tempepireset.append([])
			for i in range(len(epireset_pass[isub])):
				if len(epireset_pass[isub][i].split('|')) > 1:
					tempepireset[isub].append(epireset_pass[isub][i].split('|')[icdtime])
				else:
					tempepireset[isub].append(epireset_pass[isub][i])
		if len(gridmort_pass[isub].split('|')) > 1:
			tempgridmort.append(gridmort_pass[isub].split('|')[icdtime])
		else:
			tempgridmort.append(gridmort_pass[isub])
	
	# ----------------------
	# Variables in PopVars
	# ----------------------
	
	# Get birth rate values here, r and K_env
	if isinstance(offno, (list,tuple)):
		tempoffno = offno[icdtime]
	else:
		tempoffno = offno
	if isinstance(K_envvals, (list,tuple)):
		tempK_env = int(K_envvals[icdtime])
	else:
		tempK_env = int(K_envvals)
		
	# Check for each instance
	if isinstance(matecdmatfile, (list,tuple)):
		matecdmatfile = datadir+matecdmatfile[icdtime]
	else:	
		matecdmatfile = datadir+matecdmatfile
	if isinstance(dispcdmatfile, (list,tuple)):
		dispcdmatfile = datadir+dispcdmatfile[icdtime]
	else:
		dispcdmatfile = datadir+dispcdmatfile
	if isinstance(matemoveno, (list,tuple)):
		matemoveno = matemoveno[icdtime]		
	if isinstance(Fdispmoveno, (list,tuple)):
		Fdispmoveno = Fdispmoveno[icdtime]
	if isinstance(Mdispmoveno, (list,tuple)):
		Mdispmoveno = Mdispmoveno[icdtime]
	if isinstance(matemovethresh, (list,tuple)):
		matemovethresh = matemovethresh[icdtime]
	if isinstance(Fdispmovethresh, (list,tuple)):
		Fdispmovethresh = Fdispmovethresh[icdtime]
	if isinstance(Mdispmovethresh, (list,tuple)):
		Mdispmovethresh = Mdispmovethresh[icdtime]
	if isinstance(matemoveparA, (list,tuple)):
		matemoveparA = float(matemoveparA[icdtime])
	else:
		matemoveparA = float(matemoveparA)
	if isinstance(matemoveparB, (list,tuple)):
		matemoveparB = float(matemoveparB[icdtime])
	else:
		matemoveparB = float(matemoveparB)		
	if isinstance(matemoveparC, (list,tuple)):
		matemoveparC = float(matemoveparC[icdtime])
	else:
		matemoveparC = float(matemoveparC)
	if isinstance(FdispmoveparA, (list,tuple)):			
		FdispmoveparA = float(FdispmoveparA[icdtime])
	else:
		FdispmoveparA = float(FdispmoveparA)
	if isinstance(FdispmoveparB, (list,tuple)):
		FdispmoveparB = float(FdispmoveparB[icdtime])
	else:
		FdispmoveparB = float(FdispmoveparB)		
	if isinstance(FdispmoveparC, (list,tuple)):
		FdispmoveparC = float(FdispmoveparC[icdtime])
	else:
		FdispmoveparC = float(FdispmoveparC)		
	if isinstance(MdispmoveparA, (list,tuple)):
		MdispmoveparA = float(MdispmoveparA[icdtime])
	else:
		MdispmoveparA = float(MdispmoveparA)
	if isinstance(MdispmoveparB, (list,tuple)):
		MdispmoveparB = float(MdispmoveparB[icdtime])
	else:
		MdispmoveparB = float(MdispmoveparB)		
	if isinstance(MdispmoveparC, (list,tuple)):
		MdispmoveparC = float(MdispmoveparC[icdtime])
	else:
		MdispmoveparC = float(MdispmoveparC)
	if isinstance(twinning, (list,tuple)):
		twinning = float(twinning[icdtime])/100.
	else:
		twinning = float(twinning)/100.
	
	# Files in popvars to read in
	if isinstance(betaFile_selection,(list,tuple)):
		tempbetaFile_selection = datadir+betaFile_selection[icdtime]
	else:
		tempbetaFile_selection = datadir+betaFile_selection
	if isinstance(betaFile_epigene,(list,tuple)):
		tempbetaFile_epigene = datadir+betaFile_epigene[icdtime]
	else:
		tempbetaFile_epigene = datadir+betaFile_epigene
	
	# ---------------------------------------------
	# Read in Beta File for Multiple loci selection
	# ---------------------------------------------
	tempbetas_selection = []	
	if cdevolveans.split('_')[0] == 'M':
		# Read in Beta File
		betavals = ReadXY(tempbetaFile_selection)
	
		# Error check on beta file - should be x number of betas alleles * number of xvars * number of loci under selection
		if (len(betavals)-1)*(len(betavals[0])-1) != int(cdevolveans.split('_')[3].split('A')[1])*int(cdevolveans.split('_')[1].split('X')[1])*int(cdevolveans.split('_')[2].split('L')[1]):
			print('Beta file for selection is incorrect. Specify the beta for each variableXlociXallele combination.')
			sys.exit(-1)
		# Then extract betavals - in order of variable, loci, allele [var][loci][allele]
		for ixvar in range(int(cdevolveans.split('_')[1].split('X')[1])):
			tempbetas_selection.append([])
			for iloci in range(int(cdevolveans.split('_')[2].split('L')[1])):
				tempbetas_selection[ixvar].append([])
				for iall in range(int(cdevolveans.split('_')[3].split('A')[1])):
					tempbetas_selection[ixvar][iloci].append(float(betavals[iall+1][ixvar*(int(cdevolveans.split('_')[2].split('L')[1]))+iloci+1]))
		# Add beta0 - will be the last spot in betas vars 
		if len(betavals[0][0]) == 0:
			tempbetas_selection.append(0.0)
		else:
			tempbetas_selection.append(float(betavals[0][0])) 
	
	# ---------------------------------------------
	# Read in Beta File for Epigenetic fitness
	# ---------------------------------------------
	tempbetas_epigene = []
	if epigeneans != 'N':
		# Read in Beta File
		betavals = ReadXY(tempbetaFile_epigene)		
		# Error check on beta file - number of loci under selection * number of alleles
		if (len(betavals)-1)*(len(betavals[0])-1) != int(epigeneans.split('_')[2].split('A')[1])*int(epigeneans.split('_')[1].split('L')[1]):
			print('Beta file for epigenetics is incorrect. Specify the beta for each lociXallele combination.')
			sys.exit(-1)
		# Then extract betavals - in order of loci, allele [loci][allele]
		for iloci in range(int(epigeneans.split('_')[1].split('L')[1])):
			tempbetas_epigene.append([])
			for iall in range(int(epigeneans.split('_')[2].split('A')[1])):
				tempbetas_epigene[iloci].append(float(betavals[iall+1][iloci+1]))
		
		# Add beta0 - will be the last spot in betas vars 
		if len(betavals[0][0]) == 0:
			tempbetas_epigene.append(0.0)
		else:
			tempbetas_epigene.append(float(betavals[0][0]))	
						
	# ---------------------------------------------------------
	# Read in cdmatrix.csv and convert to a probability matrix
	# ---------------------------------------------------------
	# If mate and disp are the same, then only read in once.
	if (matecdmatfile == dispcdmatfile) and (Fdispmoveno == Mdispmoveno == matemoveno) and (Fdispmovethresh == Mdispmovethresh == matemovethresh):
		tupReadMat = ReadCDMatrix(matecdmatfile,matemoveno,\
		matemovethresh,matemoveparA,matemoveparB,matemoveparC,subpop)
		
		# Unpack tuple
		matecdmatrix = np.asarray(tupReadMat[0])
		matemovethresh = tupReadMat[1]
		mate_ScaleMin = tupReadMat[2]
		mate_ScaleMax = tupReadMat[3]
		
		# Then Set disp = mate
		Fdispcdmatrix = copy.deepcopy(matecdmatrix)
		Mdispcdmatrix = copy.deepcopy(matecdmatrix)
		Fdispmovethresh = copy.deepcopy(matemovethresh)
		Mdispmovethresh = copy.deepcopy(matemovethresh)
		Fdisp_ScaleMin = copy.deepcopy(mate_ScaleMin)
		Fdisp_ScaleMax = copy.deepcopy(mate_ScaleMax)
		Mdisp_ScaleMin = copy.deepcopy(mate_ScaleMin)
		Mdisp_ScaleMax = copy.deepcopy(mate_ScaleMax)		

	# Else if anything is different	
	else: 
		# ---------------------------------------
		# Read in cdmatrix.csv - For Mating
		# ---------------------------------------	
		tupReadMat = ReadCDMatrix(matecdmatfile,matemoveno,\
		matemovethresh,matemoveparA,matemoveparB,matemoveparC,subpop)
		matecdmatrix = np.asarray(tupReadMat[0])
		matemovethresh = tupReadMat[1]
		mate_ScaleMin = tupReadMat[2]
		mate_ScaleMax = tupReadMat[3]
	
		# ------------------------------------------------
		# Read in cdmatrix.csv - For Female Dispersal 
		# ------------------------------------------------	
		tupReadMat = ReadCDMatrix(dispcdmatfile,Fdispmoveno,\
		Fdispmovethresh,FdispmoveparA,FdispmoveparB,FdispmoveparC,subpop)
		Fdispcdmatrix = np.asarray(tupReadMat[0])
		Fdispmovethresh = tupReadMat[1]
		Fdisp_ScaleMin = tupReadMat[2]
		Fdisp_ScaleMax = tupReadMat[3]

		# ----------------------------------------------
		# Read in cdmatrix.csv - For Male Dispersal
		# ----------------------------------------------	
		tupReadMat = ReadCDMatrix(dispcdmatfile,Mdispmoveno,\
		Mdispmovethresh,MdispmoveparA,MdispmoveparB,MdispmoveparC,subpop)
		Mdispcdmatrix = np.asarray(tupReadMat[0])
		Mdispmovethresh = tupReadMat[1]
		Mdisp_ScaleMin = tupReadMat[2]
		Mdisp_ScaleMax = tupReadMat[3]
	
	# Return this functions variables
	tupClimate = matecdmatrix,Fdispcdmatrix,Mdispcdmatrix,matemovethresh,\
	Fdispmovethresh,Mdispmovethresh,Fdisp_ScaleMin,Fdisp_ScaleMax,Mdisp_ScaleMin,Mdisp_ScaleMax,mate_ScaleMin,mate_ScaleMax,tempMmort,tempFmort,tempoffno,templmbda,tempsigma,tempK_env,tempMnewmort,tempFnewmort,tempfitvals,matemoveno,Fdispmoveno,Mdispmoveno,twinning,tempMmat,tempFmat,tempbetas_selection,tempxvars_betas,tempepimod,tempepireset,tempbetas_epigene,tempgridmort	
	return tupClimate
	
	#End::DoCDClimate()
	
# ---------------------------------------------------------------------------- 
def InitializeInfect(cdinfect,Infected,nogrids,sex):
	'''
	InitializeInfect()
	This function initializes the infection status.
	'''
	
	# Create empty list to append to 
	infection = []
	tempinf = []
	
	# If cdinfect...
	if cdinfect == 'Y':
		
		# Assign initial individual a random infection status
		for ithinfect in range(nogrids):
		
			# If a NA value
			if sex[ithinfect] == 'NA':
				# Append a NA
				infection.append('NA')
			else:
				# Append a random number 0 or 1
				infection.append(int(2*np.random.uniform()))
				tempinf.append(infection[ithinfect])
	
	# If not cdinfect...
	else: 
		
		# Assign initial individual a random infection status
		for ithinfect in range(nogrids):
		
			# If a NA value
			if sex[ithinfect] == 'NA':
				# Append a NA
				infection.append('NA')
			else:			
				# Append 0
				infection.append(0)
				tempinf.append(infection[ithinfect])
			
	# Append to Infected
	Infected.append(sum(tempinf))
	del(tempinf)
	
	tupInf = infection,Infected
	return tupInf
	
	# End::InitializeInfect()

# ---------------------------------------------------------------------------- 
def DoGridOut_cdpop0(ithmcrundir,gen,loci,alleles,nogrids,subpop,xgrid,ygrid,\
id,sex,age,agelst,genes,intgenesans,infection,allelst,subpopmigration,subpopemigration,geneswap,unicor_out):
	'''
	DoGridOut_cdpop0()
	Output grid0.csv in cdpop format
	'''			
	hindex = []
	
	# Create file to write matrix to
	outputfile = open(ithmcrundir+'grid'+str(0)+'.csv','w')
	if unicor_out == 'Y' or unicor_out == True:
		outputfile_uni = open(ithmcrundir+'XY'+str(0)+'.csv','w')
		outputfile_uni.write('XCOORD,YCOORD\n')
		
	# Write out the titles 
	title = ['Subpopulation,XCOORD,YCOORD,ID,sex,age,infection,DisperseCDist,hindex,']
	outputfile.write(title[0])
	
	
	# Write out the genes title informations
	# Loop through loci
	for i in range(loci-1):
		
		# Loop for allele length
		for j in range(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	
	# To get a return character on the end of the title
	for i in range(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1))
	# Get return character
	outputfile.write('\n')
	
	# Write out all of the information.		
	for i in range(nogrids):
		outputfile.write(subpop[i]+',')
		outputfile.write(str(float(xgrid[i]))+',')
		outputfile.write(str(float(ygrid[i]))+',')		
		if sex[i] == 'NA':
			outputfile.write('OPEN,') #id
			outputfile.write('NA,') # sex
			outputfile.write('NA,') # age
			if (intgenesans != 'known'): 
				if (intgenesans != 'file_introduce') and (intgenesans != 'file_introduce_var'):
					age.append('NA') # Store age for other cases
			outputfile.write('NA,') # infection
			outputfile.write('NA,') # dispersalDist
		else:
			if unicor_out == 'Y' or unicor_out == True:
				outputfile_uni.write(str(float(xgrid[i]))+',')
				outputfile_uni.write(str(float(ygrid[i]))+'\n')
			outputfile.write(id[i]+',')
			outputfile.write(str(sex[i])+',')
			if intgenesans == 'known' or intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
				outputfile.write(age[i]+',')			
			else:
				if len(agelst) > 1: # For multiple AgeVars files
					agetemp = w_choice_general(agelst[int(subpop[i])-1])[0]
				else:
					agetemp = w_choice_general(agelst[0])[0]
				outputfile.write(str(agetemp)+',')
				age.append(agetemp)
			outputfile.write(str(infection[i])+',')
			outputfile.write('Initial,')
		
		# if known genes 
		if intgenesans == 'known':	
			if geneswap == 0:
				
				# ---------------------------------------------
				# Get AA / aa p value for genetag
				# ---------------------------------------------
				if sex[i] == 'NA':
					outputfile.write('NA,')
					hindex.append(-9999)
				else:
					if genes[i][0] == 2:
						hindex.append(1.0)
					elif genes[i][1] == 2:
						hindex.append(0.0)
					elif genes[i][0] == 1 and genes[i][1] == 1:
						hindex.append(0.5)
					else:
						hindex.append(-9999)
					outputfile.write(str(hindex[i])+',')
				
				# Write out gene info
				for iall in range(len(genes[i])):
					if sex[i] == 'NA':
						outputfile.write('NA,')
					else:
						outputfile.write(str(genes[i][iall])+',')
				outputfile.write('\n')
			else:
				# Store for hindex
				hindex.append(-9999)
							
		# These allelst values correspond to multiple subpops
		elif intgenesans == 'file' or intgenesans == 'random' or intgenesans == 'random_var' or intgenesans == 'file_var':
			if geneswap != 0: # NA values for when geneswap > gen:
				# ---------------------------------------------
				# Get AA / aa p value for genetag
				# ---------------------------------------------
				if sex[i] == 'NA':
					outputfile.write('NA,')
					hindex.append(-9999)
				else:
					if genes[i][0] == 2:
						hindex.append(1.0)
					elif genes[i][1] == 2:
						hindex.append(0.0)
					elif genes[0][0] == 1 and genes[i][1] == 1:
						hindex.append(0.5)
					else:
						hindex.append(-9999)
					outputfile.write(str(hindex[i])+',')				
				# Write out gene info
				for iall in range(len(genes[i])):
					if sex[i] == 'NA':
						outputfile.write('NA,')
					else:
						outputfile.write(str(genes[i][iall])+',')
				outputfile.write('\n')
			else:
				#Store empty array to be appended to for gene info
				indall = []							
				# And store genes information
				genes.append([])							
				# For each loci:
				for j in range(loci):
				
					# Take a random draw from the w_choice function at jth locus
					if len(allelst) > 1:
						rand1 = w_choice_general(allelst[int(subpop[i])-1][j])[0]
						rand2 = w_choice_general(allelst[int(subpop[i])-1][j])[0]
					else:					
						rand1 = w_choice_general(allelst[0][j])[0]
						rand2 = w_choice_general(allelst[0][j])[0]

					# Store genes loci spot
					#genes[i].append([])
					
					# Append assinment onto indall array - run through each condition for assignment of 1s or 2s or 0s
					# 	1s = heterozygous at that locus
					#	2s = homozygous at that locus
					#	0s = absence of allele
					for k in range(alleles[j]):
						
						# Somebody not in this spot
						if sex[i] == 'NA':
							# THen append tempindall to indall
							indall.append('NA')
							# And to genes list
							#genes[i][j].append('NA')
							genes[i].append('NA')
							
						# Else if somebody is in spot, assign genes
						else:					
							# Assignment of 2, the rest 0
							if rand1 == rand2: 
								if k < rand1 or k > rand1:
									tempindall = 0
								elif k == rand1:
									tempindall = 2
									
							# Assignment of 1s, the rest 0
							if rand1 != rand2:
								if k < min(rand1,rand2) or k > max(rand1,rand2):
									tempindall = 0
								elif k == rand1 or k == rand2:
									tempindall = 1
								else:
									tempindall = 0
									
							# THen append tempindall to indall
							indall.append(tempindall)						
							# Add to genes list
							#genes[i][j].append(tempindall)
							genes[i].append(tempindall)
				
				# ---------------------------------------------
				# Get AA / aa p value for genetag
				# ---------------------------------------------
				if sex[i] == 'NA':
					outputfile.write('NA,')
					hindex.append(-9999)
				else:
					if genes[i][0] == 2:
						hindex.append(1.0)
					elif genes[i][1] == 2:
						hindex.append(0.0)
					elif genes[i][0] == 1 and genes[i][1] == 1:
						hindex.append(0.5)
					else:
						hindex.append(-9999)
					outputfile.write(str(hindex[i])+',')
				
				# Add indall information to outputfile text
				for j in range(len(indall)-1):
					outputfile.write(str(indall[j])+',')				
				# To get return character on the end
				outputfile.write(str(indall[len(indall)-1])+'\n')
		
		# These allelst values correspond to introduced individuals through time, use first one.
		elif intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
			#Store empty array to be appended to for gene info
			indall = []							
			# And store genes information
			genes.append([])							
			# For each loci:
			for j in range(loci):
			
				# Take a random draw from the w_choice function at jth locus
				# Using the first allele file in list
				rand1 = w_choice_general(allelst[0][j])[0]
				rand2 = w_choice_general(allelst[0][j])[0]

				# Store genes loci spot
				#genes[i].append([])
				
				# Append assinment onto indall array - run through each condition for assignment of 1s or 2s or 0s
				# 	1s = heterozygous at that locus
				#	2s = homozygous at that locus
				#	0s = absence of allele
				for k in range(alleles[j]):
					
					# Somebody not in this spot
					if sex[i] == 'NA':
						# THen append tempindall to indall
						indall.append('NA')
						# And to genes list
						genes[i].append('NA')
						
					# Else if somebody is in spot, assign genes
					else:					
						# Assignment of 2, the rest 0
						if rand1 == rand2: 
							if k < rand1 or k > rand1:
								tempindall = 0
							elif k == rand1:
								tempindall = 2
								
						# Assignment of 1s, the rest 0
						if rand1 != rand2:
							if k < min(rand1,rand2) or k > max(rand1,rand2):
								tempindall = 0
							elif k == rand1 or k == rand2:
								tempindall = 1
							else:
								tempindall = 0
								
						# THen append tempindall to indall
						indall.append(tempindall)						
						# And to genes list
						genes[i].append(tempindall)
		
			# ---------------------------------------------
			# Get AA / aa p value for genetag
			# ---------------------------------------------
			if sex[i] == 'NA':
				outputfile.write('NA,')
				hindex.append(-9999)
			else:
				if genes[i][0] == 2:
					hindex.append(1.0)
				elif genes[i][1] == 2:
					hindex.append(0.0)
				elif genes[i][0] == 1 and genes[i][1] == 1:
					hindex.append(0.5)
				else:
					hindex.append(-9999)
				outputfile.write(str(hindex[i])+',')
			
			# Add indall information to outputfile text
			for j in range(len(indall)-1):
				outputfile.write(str(indall[j])+',')				
			# To get return character on the end
			outputfile.write(str(indall[len(indall)-1])+'\n')
		
	
	# Here add unique number of subpop spot in list (make list of lists)
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	# THen add spot in empty
	subpopmigration.append([])
	subpopemigration.append([])
	for i in range(nosubpops):
		subpopmigration[0].append([0])
		subpopemigration[0].append([0])
	
	# Close file
	if unicor_out == 'Y' or unicor_out == True:
		outputfile_uni.close()
	outputfile.close()
	
	# Return variables
	return genes,subpopmigration,subpopemigration,hindex
	
	# End::DoGridOut_cdpop0()	
	
# ---------------------------------------------------------------------------------------------------	 
def DoPreProcess(outdir,ibatch,ithmcrun,\
xyfilename,agefilename,equalsexratio,\
loci,intgenesans,allefreqfilename,alleles,gen,logfHndl,\
cdevolveans,cdinfect,Infected,
subpopmigration,subpopemigration,datadir,geneswap,epigeneans,unicor_out):
	'''
	DoPreProcess()
	This function does all the pre-processing work before
	CDPOP begins its time loops.
	'''
	# ----------------------------
	# Create directory
	# ----------------------------		
	ithmcrundir = outdir+'batchrun'+\
	str(ibatch)+'mcrun'+str(ithmcrun)+'/'
	os.mkdir(ithmcrundir)
	
	# ------------------------------------------------------------------
	# Read in xy points file and store info in list
	# ------------------------------------------------------------------
	xy = ReadXY(datadir+xyfilename[0]) # For PreProcessing, use the first xyfile
	
	# Error statement for 5 column data - check special cases
	# File and random use 18, file_introduce has 20, known has 21 plus genotypes, M adds more to each
	if intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
		if len(xy[1]) != 20:
			print('File_introduce option uses 20 columns in XY file, see example input files.')
			sys.exit(-1)
	elif intgenesans =='file' or intgenesans == 'random' or intgenesans == 'random_var' or intgenesans == 'file_var':
		if cdevolveans.split('_')[0] == 'M':	# Mloci selection On		
			if epigeneans == 'N': # No epigenetics
				if (len(xy[1]) - int(cdevolveans.split('_')[1].split('X')[1])) != 18:
					print('XY input file must be 18 columns plus the specified number of variables operating in the multiple loci selection model; see example input files.')
					sys.exit(-1)
			else: # epigenetics on
				if len(xy[1]) - int(cdevolveans.split('_')[1].split('X')[1]) - 2*int(epigeneans.split('_')[1].split('L')[1]) != 18:
					print('XY input file must be 17 columns plus the variables and reset parameters needed for epigenetic model, see example input files.')
					sys.exit(-1)
		elif cdevolveans.split('_')[0] == 'Hindex': #Hindex option
			if epigeneans != 'N':
				print('Epigenetics and Hindex options currently not implemented together.')
				sys.exit(-1)
			else:
				if len(xy[1]) != 19: # check the length on these
					print('Hindex option given for cdevolveans and XY input file must be 19 columns; see example input files.')
					sys.exit(-1)
				
		else: # anything else
			if 	epigeneans == 'N': # No epigenetics
				if len(xy[1]) != 18:
					print('File and random options use 18 columns in XY file, see example input files.')
					sys.exit(-1)
			else: # epigenetics
				if len(xy[1]) - 2*int(epigeneans.split('_')[1].split('L')[1]) != 18:
					print('XY input file must be 18 columns plus the variables and reset parameters needed for epigenetic model; see example input files.')
					sys.exit(-1)
	
	# Store all information in lists by variable name
	FID = []
	subpop = []
	xgrid = []
	ygrid=[]
	id = []
	fitvals = [] # selection for 1-2 locus model
	xvars = [] # selection for multilocus/hindex options 
	sex = []	
	age = []	
	genes = []
	infection = []
	epimod_vals = [] # epigenetic modification
	epireset_vals = [] # epigenetic reset values
	hindex = [] # Hindex value
	gridmort = [] # Added mortality for each grid location
	
	for i in range(len(xy)-1):
		FID.append(i)
		subpop.append(xy[i+1][0])
		if xy[i+1][1] == 'NA' or xy[i+1][2] == 'NA':
			print('You need to specify the (x,y) locations for this location even if it is a NA value.')
			sys.exit(-1)
		if xy[i+1][3] == 'NA' and equalsexratio == 'WrightFisher':
			print('You can not force equal sex ratio (WrightFisher) with fluctuating population size.')
			print('Specify N for Equal Sex Ratio option.')
			sys.exit(-1)
		xgrid.append(float(xy[i+1][1]))
		ygrid.append(float(xy[i+1][2]))
		gridmort.append(xy[i+1][3])
		id.append(xy[i+1][4])
		
		# Only grab sex information from the file if equal sex ratio is N
		# -------------------------
		if equalsexratio == 'N' or equalsexratio == 'AtBirth':
			if xy[i+1][3] == 'NA' and xy[i+1][5] != 'NA':
				print('Must specify NA for the sex value of all empty starting locations.')
				sys.exit(-1)
			sex.append(xy[i+1][5])
				
			# Change id to 'OPEN'
			if sex[i] == 'NA':
				id[i] = 'OPEN'
			
			# If sex was F or M...change to 0 and 1
			if sex[i] == 'F':
				sex[i] = '0'
			elif sex[i] == 'M':
				sex[i] = '1'
		# ---------------------------
		
		# Get fitvals and xvars - indexspot is the start of fitvals
		# ---------------------------------------------------------
		# For file or random case - get index spots for fitness values (1-2 locus model)
		if intgenesans == 'file' or intgenesans == 'random' or intgenesans == 'file_var' or intgenesans == 'random_var': # 17 columns 
			indexspot = 6
		# for known case
		elif intgenesans == 'known':
			indexspot = 9
		# for file_introduce
		elif intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
			indexspot = 8
		else:
			print('intgenesans does not exist, see user manual for options.')
			sys.exit(-1)
		xvars_indexspot = indexspot + 12 # for getting the xvars for M and Hindex options - fitness spots are length 12 always.
		
		# Get fitvals and xvars here
		if cdevolveans == '1' or cdevolveans == '3' or cdevolveans == '1_HeMort_GEA' or cdevolveans == '1_HeMort_All':
			fitvals.append([xy[i+1][indexspot+0],xy[i+1][indexspot+1],xy[i+1][indexspot+2]])
		elif cdevolveans == '2':
			fitvals.append([xy[i+1][indexspot+3],xy[i+1][indexspot+4],xy[i+1][indexspot+5],xy[i+1][indexspot+6],xy[i+1][indexspot+7],xy[i+1][indexspot+8],xy[i+1][indexspot+9],xy[i+1][indexspot+10],xy[i+1][indexspot+11]])
		elif cdevolveans.split('_')[0] == 'M':
			xvars.append([])
			for ixvars in range(int(cdevolveans.split('_')[1].split('X')[1])):
				xvars[i].append(xy[i+1][xvars_indexspot+ixvars])
		elif cdevolveans.split('_')[0] == 'Hindex':
			xvars.append([])# So far just 1 xvars can be given for Hindex option
			xvars[i].append(xy[i+1][xvars_indexspot])
				
		# Get epigenetic information, comes after selection, 
		if epigeneans != 'N':
			epimod_vals.append([])
			epireset_vals.append([])
			if cdevolveans.split('_')[0] != 'M':
				for ixvars in range(int(epigeneans.split('_')[0].split('X')[1])):
					epimod_vals[i].append(xy[i+1][xvars_indexspot+ixvars])
				for ireset in range(int(epigeneans.split('_')[1].split('L')[1])):
					epireset_vals[i].append(xy[i+1][xvars_indexspot+int(epigeneans.split('_')[0].split('X')[1])+ireset])
			else:
				for ixvars in range(int(epigeneans.split('_')[0].split('X')[1])):
					epimod_vals[i].append(xy[i+1][xvars_indexspot+ixvars+int(cdevolveans.split('_')[1].split('X')[1])])
				for ireset in range(int(epigeneans.split('_')[1].split('L')[1])):
					epireset_vals[i].append(xy[i+1][xvars_indexspot+int(epigeneans.split('_')[0].split('X')[1])+ireset+int(cdevolveans.split('_')[1].split('X')[1])])
			'''
			# Leaving this indexing code here when I had epistasis region
			elif cdevolveans.split('_')[0] == 'M' and epistasis != 'N':
				for ixvars in xrange(int(epigeneans.split('_')[0].split('X')[1])):
					epimod_vals[i].append(xy[i+1][xvars_indexspot+1+ixvars+int(cdevolveans.split('_')[1].split('X')[1])])
				for ireset in xrange(int(epigeneans.split('_')[1].split('L')[1])):
					epireset_vals[i].append(xy[i+1][xvars_indexspot+1+int(epigeneans.split('_')[0].split('X')[1])+ireset+int(cdevolveans.split('_')[1].split('X')[1])])
			'''
		# -----------------------------------------------------------------------------------------		
		
		# Grab age and genes for special cases
		# ----------------------------------------------------------
		if intgenesans == 'known':
			knownindexspot = 22 # index spot for start of genes
			if geneswap != 0:# switch to random genes			
				print('Warning: Known gene file can not be used with gene start time > 0 Will continue but with random gene assignment.')
				intgenesans = 'random'
			else:
				addindex = 0 # in case selection or epigenetics added
				if epigeneans != 'N':
					addindex = addindex + int(epigeneans.split('_')[0].split('X')[1]) + int(epigeneans.split('_')[1].split('L')[1])
				if cdevolveans.split('_')[0] == 'M':
					addindex = addindex + int(cdevolveans.split('_')[1].split('X')[1])
				
				if len(xy[1]) != knownindexspot+addindex+sum(alleles):
					print('Specified known genetic initializtion file. Not correct format. See example files for number of fields needed.')
					sys.exit(-1)
				else:					
					# Age storage here for known file
					age.append(xy[i+1][6])
					# Error check for age in known file must be 1 or greater
					if age[i] != 'NA':
						if int(age[i]) < 1:
							print('Known file must initize with age 1+.')
							sys.exit(-1)
					# genes[individual]
					#genes.append([])			
					# Error check here to make sure gene file matches specified loci and alleles
					if sum(alleles) != len(xy[i+1][int(knownindexspot+addindex):len(xy[i+1])]):
						print('Known genes file does not match loci and alleles given.')
						sys.exit(-1)							
					genes.append(xy[i+1][knownindexspot+addindex:])
			
		# For other cases when geneswap is great than 0
		if intgenesans != 'known' and geneswap != 0:
			genes.append(['NA']*(sum(alleles))) # genes[individual]
			
		# -----------------------------------------------
		# Grab age and infection for 'file_introduce' special case
		if intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
			# Age storage here for known file
			age.append(xy[i+1][6])
			# Error check for age in known file must be 1 or greater
			if age[i] != 'NA':
				if int(age[i]) < 1:
					print('Known file must initize with age 1+.')
					sys.exit(-1)
			infection.append(xy[i+1][7])
	
	# Tally infection
	if intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':		
		Infected.append(len(np.where(np.asarray(infection)=='1')[0]))		
	# Store the number of grids
	nogrids = len(xy)-1
	
	# Delete x variable
	del(xy)
		
	# --------------------------
	# Error Checks
	# --------------------------	
	# For now, subpops need to be ordered 1 to N and not skipping, no 0s
	if len(np.where(np.unique(subpop)=='0')[0]) != 0:
		print('Subpopulation identification field can not have 0 values.')
		sys.exit(-1)
	tempcheck = []
	for i in range(len(np.unique(subpop))):
		tempcheck.append(int(np.unique(subpop)[i]))
	tempcheck = np.sort(tempcheck)
	if len(tempcheck) > 1:
		for i in range(len(tempcheck)-1):
			if tempcheck[i+1]-tempcheck[i] > 1:
				print('Subpopulation identification field must be labeled sequentially or a single value.')
				sys.exit(-1)	
		
	# If equal sex ratio is Y, then split up sex into equal parts
	if equalsexratio == 'WrightFisher':
	
		# Get unique number of subpops, make equal sex ratio for each subpopulation
		nosubpops = len(np.unique(subpop))
		unique_subpops = np.unique(subpop)
		
		# Count up the unique number of subgrids appending to subgrids
		subgridtotal = []
		# Create list of lists storage spots for number of subgrids
		for i in range(nosubpops):
			subgridtotal.append([])
		for i in range(len(subpop)):
			# Loop through unique subpops
			for j in range(nosubpops):
				# If subpop exits append to subgrid spot
				if subpop[i] == unique_subpops[j]:
					subgridtotal[int(unique_subpops[j])-1].append(1)
		
		# And then sum them up
		for i in range(nosubpops):
			subgridtotal[i] = sum(subgridtotal[i])
			# If the subpopulation number is not even then sys exit
			if np.mod(subgridtotal[i],2) == 1:
				print("You have equal sex ratio turned on and this population is not even.")
				sys.exit(-1)
			
			# Then create half males and females and shuffle
			sextemp = np.append(np.zeros(subgridtotal[i]/2,"int"),np.ones(subgridtotal[i]/2,"int"))
			np.random.shuffle(sextemp)
			# Loop through these individuals and append to LIST sex
			for j in range(len(sextemp)):			
				# The add them together and shuffle and append to sex
				sex.append(str(sextemp[j]))
				
		# Delete extra stuff
		del(sextemp)
		
	# -------------------------------------------
	# Initialize age structure
	# ------------------------------------------- 
	agelst,Magemort,Fagemort,egg_lmbdavals,egg_sigmavals,Mnewmortperc,Fnewmortperc,Mmature,Fmature = InitializeAge(agefilename,nogrids,datadir)
	
	# ------------------------------------------
	# Initialize infection
	# ------------------------------------------
	if intgenesans != 'file_introduce' or intgenesans == 'file_introduce_var':
		tupInf = InitializeInfect(cdinfect,Infected,nogrids,sex)
		infection = tupInf[0]
		Infected = tupInf[1]
	
	# ----------------------------------------------
	# Initialize genetic structure
	# ----------------------------------------------
	allelst = InitializeGenes(intgenesans,allefreqfilename,loci,alleles,datadir)
	
	# --------------------------------------------------------------------
	# Create output file grid0.csv and write to it and return genes
	# -------------------------------------------------------------------- 
	genes,subpopmigration,subpopemigration, hindex = DoGridOut_cdpop0(ithmcrundir,0,loci,alleles,\
	nogrids,subpop,xgrid,ygrid,id,sex,age,agelst,genes,intgenesans,\
	infection,allelst,subpopmigration,subpopemigration,geneswap,unicor_out)
			
	# Return this functions variables
	tupPreProcess = ithmcrundir,FID,id,sex,age,xgrid,ygrid,genes,\
	nogrids,subpop,fitvals,infection,Infected,subpopmigration,subpopemigration,Magemort,Fagemort,egg_lmbdavals,egg_sigmavals,allelst,Mnewmortperc,Fnewmortperc,Mmature,Fmature,intgenesans,xvars,epimod_vals,epireset_vals,hindex, gridmort
	return tupPreProcess
	
	#End::DoPreProcess()
	
# ---------------------------------------------------------------------------------------------------	 		
def DoUserInput(fileans):
	
	'''
	DoUserInput()
	This function reads in the user input and 
	stores the variables.
	'''
	
	# Open file for reading
	inputfile = open(fileans,'r')

	# Read lines from the file
	lines = inputfile.readlines()

	#Close the file
	inputfile.close()

	# Create an empty matrix to append to
	inputvariables = []

	# Split up each line in file and append to empty matrix, x
	for i in lines:
		thisline = i.split(',')
		inputvariables.append(thisline)
		
	# Delete lines
	del(lines)

	return inputvariables
	
	#End::DoUserInput()