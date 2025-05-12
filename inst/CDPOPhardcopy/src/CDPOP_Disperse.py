# -------------------------------------------------------------------------------------------------
# CDPOP_Disperse.py
# Author: Erin L Landguth
# Created: December 2010
# Description: This is the function/module file for dispersal processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError("Numpy required.")

# Python specific functions
import pdb, random, copy, os, sys

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False

# --------------------------------------------------------------------------
def logMsg(outf,msg):
	'''
	logMsg() --log file message handler.
	Inputs:
	outf - open file handle
	msg -- string containing formatted message
	--always outputs to log file by default.
	--using msgVerbose, can be set to "Tee" output to stdout as well
	'''
	outf.write(msg+ '\n')
	if msgVerbose:
		print(("%s"%(msg)))
		
	# End::logMsg()

# ---------------------------------------------------------------------------------------------------	 
def w_choice_general(lst):
	'''
	w_choice_general()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(x[1] for x in lst)
	n=np.random.uniform(0,wtotal)
	count = 0
	for item, weight in lst:
		if n < weight:
			break
		n = n-weight
		count = count + 1
	# The case where all of the values in lst are the same
	if len(lst) == count:
		count = count-1
	return item,count
	
	#End::w_choice_general()	

# ---------------------------------------------------------------------------------------------------	 
def w_choice_item(lst):
	'''
	w_choice_item()
	Weighted random draw from a list, probilities do not have to add to one.
	'''
	wtotal=sum(lst)
	n=np.random.uniform(0,wtotal)
	for i in range(len(lst)):
		if n < lst[i]:
			break
		n = n-lst[i]
	return i
	
	#End::w_choice_item()
	
# ---------------------------------------------------------------------------------------------------	
def GetProbArray(Fxycdmatrix,Mxycdmatrix,tempoffspring,index,freegrid,philopatry,F_freegrid=None,M_freegrid=None,F_off=None,M_off=None,Fcount=None,Mcount=None):
	'''
	GetProbArray()
	This function gets indices for F and M specific cdmatrix values
	'''
		
	# No philopatry
	if philopatry == 'N':
	
		# Index into offspring array
		currentoff = tempoffspring[index]
		
		# Append the freegrid probabilities for the offspring choice
		if int(currentoff[4]) == 0: # Female offspring
			probarray = Fxycdmatrix[int(currentoff[0])][freegrid]
			Fcount = Fcount + 1
			sexans = 'F'
		elif int(currentoff[4]) == 1: # Male offspring
			probarray = Mxycdmatrix[int(currentoff[0])][freegrid]
			Mcount = Mcount + 1
			sexans = 'M'
		else:
			print('Invalid offspring list.')
			sys.exit(-1)	

	elif philopatry == 'F':
		
		# Get the first half of freegrid filled up by females
		if Fcount < min(F_freegrid,len(F_off)):
			probarray = Fxycdmatrix[int(F_off[Fcount][0])][freegrid]
			Fcount = Fcount + 1
			sexans = 'F'
		elif Mcount < min(M_freegrid,len(M_off)):
			probarray = Mxycdmatrix[int(M_off[Mcount][0])][freegrid]
			Mcount = Mcount + 1
			sexans = 'M'
		# Extra spots and offspring
		else:
			# Randomly select the remaining
			females_left = len(F_off) - Fcount
			if females_left <= 0:
				females_left = np.zeros(0,dtype = int)
			else:
				females_left = np.zeros(females_left,dtype = int)
			males_left = len(M_off) - Mcount
			if males_left <= 0:
				males_left = np.ones(0,dtype = int)
			else:
				males_left = np.ones(males_left,dtype = int)
			all_left_sex = np.concatenate([females_left,males_left])
			
			if len(all_left_sex) != 0:
				rand_left = np.random.randint(0, len(all_left_sex) - 1)
				if all_left_sex[rand_left] == 0:
					probarray = Fxycdmatrix[int(F_off[Fcount][0])][freegrid]
					Fcount = Fcount + 1
					sexans = 'F'
				else:
					probarray = Mxycdmatrix[int(M_off[Mcount][0])][freegrid]
					Mcount = Mcount + 1
					sexans = 'M'
			else:
				probarray = np.zeros(0,dtype=float)
				sexans = 'N'
				
		
	elif philopatry == 'M':		
		
		# Get the first half of freegrid filled up by males
		if Mcount < min(M_freegrid,len(M_off)):
			probarray = Mxycdmatrix[int(M_off[Mcount][0])][freegrid]
			Mcount = Mcount + 1
			sexans = 'M'
		elif Fcount < min(F_freegrid,len(F_off)):
			probarray = Fxycdmatrix[int(F_off[Fcount][0])][freegrid]
			Fcount = Fcount + 1
			sexans = 'F'
		# Extra spots and offspring -randomly choose the left over offspring
		else:
			# Randomly select the remaining
			females_left = len(F_off) - Fcount
			if females_left <= 0:
				females_left = np.zeros(0,dtype = int)
			else:
				females_left = np.zeros(females_left,dtype = int)
			males_left = len(M_off) - Mcount
			if males_left <= 0:
				males_left = np.ones(0,dtype = int)
			else:
				males_left = np.ones(males_left,dtype = int)
			all_left_sex = np.concatenate([females_left,males_left])
			if len(all_left_sex) != 0:	
				rand_left = np.random.randint(0, len(all_left_sex) - 1)
				if all_left_sex[rand_left] == 0:
					probarray = Fxycdmatrix[int(F_off[Fcount][0])][freegrid]
					Fcount = Fcount + 1
					sexans = 'F'
				else:
					probarray = Mxycdmatrix[int(M_off[Mcount][0])][freegrid]
					Mcount = Mcount + 1
					sexans = 'M'
			else:
				probarray = np.zeros(0,dtype=float)
				sexans = 'N'
	
	return probarray,Fcount,Mcount,sexans
	
	# End::GetProbArray()

# ---------------------------------------------------------------------------------------------------	
def DoHindexSelection(offspring,tempfreegrid,iteminlist,cdevolveans,xvars):
	'''
	This function calculates offspring differential mortality, given the individuals hindex, and pars given
	'''
	
	Hindex = offspring[-1] # hindex of offspring
	
	if Hindex != -9999:
		# Linear
		# ------	
		if cdevolveans.split('_')[1] == 'Linear':
			# Get parameters
			pars = cdevolveans.split('_')[2].split(';')
			slope_min = float(pars[0])
			slope_max = float(pars[1])
			int_min = float(pars[2])
			int_max = float(pars[3])
			X_min = float(pars[4])
			X_max = float(pars[5])
			X_val = float(xvars[tempfreegrid[iteminlist]][0])
						
			# Calculate slope
			m = ((slope_min - slope_max) / (X_min - X_max)) * X_val - X_min * ((slope_min - slope_max) / (X_min - X_max)) + slope_min
			
			# Calculate intercept
			b = ((int_max - int_min) / (X_min - X_max)) * X_val - X_min * ((int_max - int_min) / (X_min - X_max)) + int_max
			
			# Get fitness value
			Fitness = m * Hindex + b					
		
		# Error
		# -----
		else:
			print('Hindex cdevolve answer not given.')
			sys.exit(-1)
	
	# Convert fitness to differential mortality - 1 - Finess
	differentialmortality = 1. - Fitness
		
	return differentialmortality	
	# End::DoHindexSelection
	
# ---------------------------------------------------------------------------------------------------	
def DoMLocusSelection(offspring,tempfreegrid,iteminlist,loci,cdevolveans,betas,xvars,maxfit,minfit,gen):
	'''
	DoMLocusSelection()
	This function calculates offsprings differential mortality, given the individuals genotype, betas, and xvariables supplied for the linear additive model. 
	'''
	#pdb.set_trace() # Check var alleles here too?
	# Get individuals genes under selection
	# -----------------------------
	offgenes = np.asarray(offspring[7]) # Genes
	#offgenes = offgenes.tolist()
	#offgenes = [offgenes[x:x+(len(offgenes)/loci)] for x in xrange(0, len(offgenes), (len(offgenes)/loci))]
	# Assume the first N loci are under selection (and the first A alleles under selection)
	indexto = int(cdevolveans.split('_')[2].split('L')[1]) * int(cdevolveans.split('_')[3].split('A')[1])
	selgenes = offgenes[0:indexto]
	
	# Apply linear model
	atspot_xvars = xvars[tempfreegrid[iteminlist]] # X vars at this grid spot
	
	linmodel = [] # [xvar][loci][allele]
	# Loop through each environment
	for ixvar in range(int(cdevolveans.split('_')[1].split('X')[1])):
		xvar = float(atspot_xvars[ixvar])
		thesebetas = sum(betas[ixvar],[])
		# Loop through the alleles
		for iall in range(len(selgenes)):
			thisbeta = thesebetas[iall]
			thisallele = selgenes[iall]
			if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
				if thisallele == 2:
					thisallele = 1 # Change second copy to 1
			# allele X environment X beta effect
			linmodel.append(xvar * thisbeta * thisallele)		
	
	# Add the beta not
	linmodel.append(betas[-1])
	
	'''
	# For max/min options, get global vs local options
	if cdevolveans.split('_')[5] == 'Global':
		#pdb.set_trace() # maxfit[-1] check this location index
		thismaxfit = maxfit[0] 
		thisminfit = minfit[0]
	else:
		thismaxfit = maxfit[gen]
		thisminfit = minfit[gen]
	'''		
	# Get the max and min values for rescaling 
	thismaxfit = maxfit[0] 
	thisminfit = minfit[0]
	
	# For the case in which the population fixated and all the same genotypes
	if thismaxfit - thisminfit == 0:
		Fitness = 1.
	else:
		# Rescale GXE calculation by max/min 
		Fitness = (sum(linmodel) - thisminfit) / (thismaxfit - thisminfit)
	
	# Fitness could be less than 0 because of rescaling cases sum(linmodel) < thisminfit
	if Fitness < 0:
		Fitness = 0.
	# Fitness could be greater than 1 because of rescale cases sum(linmodel) > thismaxfit
	if Fitness > 1:
		Fitness = 1.
	
	# Add the linear model together and logit
	#Fitness = np.exp(sum(linmodel)) / (1. + np.exp(sum(linmodel)))
	
	# Convert fitness to differential mortality - 1 - Finess
	differentialmortality = 1. - Fitness
		
	return differentialmortality
	
	# End::DoMLocusSelection()


# ---------------------------------------------------------------------------------------------------	
def DoHeMortSelection(offspring,fitvals1,tempfreegrid,iteminlist,loci,cdevolveans):
	'''
	DoHeMortSelection()
	This function calculates offsprings differential mortality, given the individuals Het and equation specified. 
	'''
	
	# Get individuals heterozygosity - # het loci / loci
	# -----------------------------
	offgenes = np.asarray(offspring[7],dtype=int) # Genes
	#all_freq_sq = (offgenes / 2.)**2 # Allele frequency ^ 2
	#homozygosity = sum(all_freq_sq)/loci
	#he = (1. - homozygosity)	
	hetLoci = loci - len(np.where(offgenes == 2)[0])
	he = hetLoci / float(loci)
	
	# For GEA method
	if cdevolveans == '1_HeMort_GEA':	
		# If L0A0|L0A0 -- loci under selection:
		if offgenes[0] == 2:

			# The grab it's fitness values
			m = float(fitvals1[tempfreegrid[iteminlist]][0][0])
			bint = float(fitvals1[tempfreegrid[iteminlist]][0][1])
			y = m * he + bint # This is survival
			if y > 100:
				y = 100.
			differentialmortality = (100. - y)/100.
																	
		# If L0A0|L0A1 -- loci under selection:
		elif offgenes[0] == 1 and offgenes[1] == 1:

			# The grab it's fitness values
			m = float(fitvals1[tempfreegrid[iteminlist]][1][0])
			bint = float(fitvals1[tempfreegrid[iteminlist]][1][1])
			y = m * he + bint # This is survival
			if y > 100:
				y = 100.
			differentialmortality = (100. - y)/100.
														# If L0A1|L0A1 -- loci under selection
		elif offgenes[1] == 2:
			
			# The grab it's fitness values
			m = float(fitvals1[tempfreegrid[iteminlist]][2][0])
			bint = float(fitvals1[tempfreegrid[iteminlist]][2][1])
			y = m * he + bint # This is survival
			if y > 100:
				y = 100.
			differentialmortality = (100. - y)/100.
		
		# Another genotype
		else:
			differentialmortality = 0.0
		
	# For all method
	else:
		# The grab it's fitness values
		m = float(fitvals1[tempfreegrid[iteminlist]][0][0])
		bint = float(fitvals1[tempfreegrid[iteminlist]][0][1])
		y = m * he + bint # This is survival
		if y > 100:
			y = 100.
		differentialmortality = (100. - y)/100.
	
	return differentialmortality
	
	# End::DoHeMortSelection()
	
	
# ---------------------------------------------------------------------------------------------------	
def Do1LocusSelection(offspring,fitvals1,tempfreegrid,iteminlist):
	'''
	Do1LocusSelection()
	This function calculates offsprings differential mortality, ie,
	offspring viability selection, for the 1-locus selection model.
	'''
	
	# If L0A0|L0A0 -- loci under selection:
	if int(offspring[7][0]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals1[tempfreegrid[iteminlist]][0])/100.
																
	# If L0A0|L0A1 -- loci under selection:
	elif int(offspring[7][0]) == 1 and int(offspring[7][1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals1[tempfreegrid[iteminlist]][1])/100.
																															
	# If L0A1|L0A1 -- loci under selection
	elif int(offspring[7][1]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals1[tempfreegrid[iteminlist]][2])/100.
	
	# Another genotype
	else:
		differentialmortality = 0.0
	
	return differentialmortality
	
	# End::Do1LocusSelection()
	
# ---------------------------------------------------------------------------------------------------	
def Do2LocusSelection(offspring,fitvals2,tempfreegrid,iteminlist,alleles):
	'''
	Do2LocusSelection()
	This function calculates offsprings differential mortality, ie,
	offspring viability selection, for the 2-locus selection model.
	'''
	
	# If L0A0|L0A0|L1A0|L1A0 - AABB -- loci under selection:
	if int(offspring[7][0]) == 2 and int(offspring[7][alleles[0]]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][0])/100.
									
	# If L0A0|L0A1|L1A0|L1A0 - AaBB -- loci under selection:
	elif int(offspring[7][0]) == 1 and int(offspring[7][1]) == 1 and int(offspring[7][alleles[0]]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][1])/100.
																															
	# If L0A1|L0A1|L1A0|L1A0 - aaBB -- loci under selection
	elif int(offspring[7][1]) == 2 and int(offspring[7][alleles[0]]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][2])/100.
									
	# If L0A0|L0A0|L1A0|L1A1 - AABb -- loci under selection:
	elif int(offspring[7][0]) == 2 and int(offspring[7][alleles[0]]) == 1 and int(offspring[7][alleles[0]+1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][3])/100.
									
	# If L0A0|L0A1|L1A0|L1A1 - AaBb -- loci under selection:
	elif int(offspring[7][0]) == 1 and int(offspring[7][1]) == 1 and int(offspring[7][alleles[0]]) == 1 and int(offspring[7][alleles[0]+1]) == 1:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][4])/100.
																	
	# If L0A1|L0A1|L1A0|L1A1 - aaBb -- loci under selection
	elif int(offspring[7][1]) == 2 and int(offspring[7][alleles[0]]) == 1 and int(offspring[7][alleles[0]+1]) == 1:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][5])/100.
	
	# If L0A0|L0A0|L1A1|L1A1 - AAbb -- loci under selection:
	elif int(offspring[7][0]) == 2 and int(offspring[7][alleles[0]+1]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][6])/100.
									
	# If L0A0|L0A1|L1A1|L1A1 - Aabb -- loci under selection:
	elif int(offspring[7][0]) == 1 and int(offspring[7][1]) == 1 and int(offspring[7][alleles[0]+1]) == 2:

		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][7])/100.																													
	# If L0A1|L0A1|L1A1|L1A1 - aabb -- loci under selection
	elif int(offspring[7][1]) == 2 and int(offspring[7][alleles[0]+1]) == 2:
		
		# The grab it's fitness values
		differentialmortality = float(fitvals2[tempfreegrid[iteminlist]][8])/100.
	
	# Another genotype
	else:
		differentialmortality = 0.0	
		
	return differentialmortality
	
	# End::Do2LocusSelection()

# ---------------------------------------------------------------------------------------------------	
def DoEmigration(offspring,freegrid,Migrants,Open,loci,alleles,\
Fxycdmatrix,Mxycdmatrix,gen,\
offspringno,cdevolveans,fitvals,subpop,subpopmigration,DisperseDeaths,CouldNotDisperse,\
gridmort,philopatry,females,subpopemigration,females_nomate,\
males,males_nomate,startSelection,betas,xvars_betas,maxfit,minfit):
	'''
	DoEmigration()
	This function enforces emigration when there are
	more offspring than open grid spots.
	'''	
	
	# Create variable to store offspring that disperse inside grid
	OffDisperseIN=[]
	
	# Do a shuffle on offpspring
	shuffle(offspring)
		
	# Initialize the while loop
	dispcount = 0
	offcount = 0
	Fcount = 0 # Only used for philopatry options
	Mcount = 0 # Only used for philopatry options
	F_freegrid = len(freegrid)/2 # Only used for philopatry options
	M_freegrid = len(freegrid)/2 # Only used for philopatry options
	F_off = [] # Only used for philopatry options
	M_off = [] # Only used for philopatry options
	
	# Add spot to track dispersing deaths for cdevolve
	DisperseDeaths.append([])
	CouldNotDisperse.append([])
		
	# Deep copy the freegrid spots to delete from
	tempfreegrid = copy.deepcopy(freegrid)
	#pdb.set_trace()
	# If philopatry, order the offspring grid by F/M or M/F
	if philopatry == 'F' or philopatry == 'f' or philopatry == 'female' or philopatry == 'Female':
		offspring.sort(key=lambda x: x[4],reverse=False)
		# Get offspring sex
		F_off = []
		M_off = []
		for row in offspring:
			# Female offspring list
			if int(row[4]) == 0:
				F_off.append(row)
			# Male offspring list
			else:
				M_off.append(row)
						
	elif philopatry == 'M' or philopatry == 'm' or philopatry == 'Male' or philopatry == 'male':
		offspring.sort(key=lambda x: x[4],reverse=True)
		# Get offspring sex
		for row in offspring:
			# Female offspring list
			if int(row[4]) == 0:
				F_off.append(row)
			# Male offspring list
			else:
				M_off.append(row)	
	
	# This while loop makes sure loop stops at carrying capacity (ie, total number of freegrids)
	#	or stops at end of offpsring list
	while dispcount < len(freegrid) and offcount < len(offspring):
		
		# Loop through offspring that are shuffled
		for i in range(len(offspring)):
			#pdb.set_trace()	
			# Create a function here that gets indices for male and female
			probarray,Fcount,Mcount,sexans = GetProbArray(Fxycdmatrix,Mxycdmatrix,offspring,i,\
			tempfreegrid,philopatry,F_freegrid,M_freegrid,F_off,M_off,Fcount,Mcount)
								
			# If statement to check if there are spots for offpsring to disperse to
			if sum(probarray) != 0.0:
				
				# Select the w_choice item
				iteminlist = w_choice_item(probarray)
				
				# Then Check and Calculated Differential Mortality Options:
				# ---------------------------------------------------------
				
				# CDEVOLVE - 1 loci
				if cdevolveans == '1' and gen >= startSelection:
					
					# Call 1-locus selection model
					differentialmortality = Do1LocusSelection(offspring[i],fitvals,tempfreegrid,iteminlist)																		
				# CDEVOLVE - 2 loci
				elif cdevolveans == '2' and gen >= startSelection:
					
					# Call 2-locus selection model
					differentialmortality = Do2LocusSelection(offspring[i],fitvals,tempfreegrid,iteminlist,alleles)			
						
				# CDEVOLVE - Heterozygosity survival
				elif (cdevolveans == '1_HeMort_GEA' or cdevolveans == '1_HeMort_All') and gen >= startSelection:
					
					# Call 1-locus selection model
					differentialmortality = DoHeMortSelection(offspring[i],fitvals,tempfreegrid,iteminlist,loci,cdevolveans)
				
				# CDEVOLVE - Multiple loci selection model
				elif cdevolveans.split('_')[0] == 'M' and gen >= startSelection:
					
					# Call M-locus selection model
					differentialmortality = DoMLocusSelection(offspring[i],tempfreegrid,iteminlist,loci,cdevolveans,betas,xvars_betas,maxfit,minfit,gen)	

				# CDEVOLVE - Hindex
				elif (cdevolveans.split('_')[0] == 'Hindex') and (gen >= startSelection):
					 
					# Call Hindex selection model
					differentialmortality = DoHindexSelection(offspring[i],tempfreegrid, iteminlist, cdevolveans,xvars_betas)
					
				# CDEVOLVE - not on
				else:
					
					differentialmortality = 0.0
				
				# The Check for spatial mortality - subpopulation mortality				
				#fromsubpop = subpop[int(offspring[i][0])] # What subpopulation is offspring coming from
				tosubpop = subpop[tempfreegrid[iteminlist]] # Where is subpopulation proposing to go
				
				# Grab its mortality percentage
				#differentialmortality_SpatialSubPopMort = float(gridmort[int(tosubpop)-1][0])/100.
				differentialmortality_SpatialSubPopMort = float(gridmort[tempfreegrid[iteminlist]])/100.
				
				# Calculated and Check Differential Mortality
				# -------------------------------------------
				differentialmortality_Total = 1. - ((1. - differentialmortality) * (1. - differentialmortality_SpatialSubPopMort))
				
				# Then flip the coin to see if offspring survives its location
				randcheck = np.random.uniform()
				
				# If offspring did not survive: break from loop, move to next offspring
				if randcheck < differentialmortality_Total:
					offcount = offcount + 1
					DisperseDeaths[gen].append(1)
					CouldNotDisperse[gen].append(0)
					continue
								
				# Append information to variable [offspring, grid it dispersed to, and name]
				recd = [offspring[i],tempfreegrid[iteminlist],'T'+str(gen)+\
				'M'+str(offspring[i][0])+'F'+str(offspring[i][1])+\
				'Pop'+str(subpop[int(offspring[i][0])])]
							
				# Record offspring disperse information	
				OffDisperseIN.append(recd)
									
				# Update count for freegrid filling up
				dispcount = dispcount + 1
				offcount = offcount + 1
				DisperseDeaths[gen].append(0)
				CouldNotDisperse[gen].append(0)
								
				# Store the subpop dispersing to another subpop...what if 1, 3, 5 labeled subpop. either error intially or fix here....
				dispersingfrom = subpop[int(recd[0][0])]
				dispersingto = subpop[tempfreegrid[iteminlist]]
				if dispersingto != dispersingfrom:
					subpopmigration[gen][int(dispersingto)-1].append(1)
					subpopemigration[gen][int(dispersingfrom)-1].append(1)
					
				# And then delete freegrid spot from temp variable
				del(tempfreegrid[iteminlist])
				
			# If statement to check if there were not spots to disperse to
			elif sum(probarray) == 0.0:
				
				# Then Break from the loop and move to next offspring
				offcount = offcount + 1
				CouldNotDisperse[gen].append(1)
				
				# For philopatry options and tracking only
				if philopatry != 'N':
					if sexans == 'F':
						Fcount = Fcount + 1
					else:
						Mcount = Mcount + 1
				
				continue
										
			
	# The store the number of disperses to separate subpops
	subpopmigration.append([]) # This adds a spot for next generation
	subpopemigration.append([])
	for i in range(len(subpopmigration[0])):
		subpopmigration[gen][i]=sum(subpopmigration[gen][i])
		subpopmigration[gen+1].append([0]) # This adds spots for subpops in next generation
		subpopemigration[gen][i]=sum(subpopemigration[gen][i])
		subpopemigration[gen+1].append([0]) # This adds spots for subpops in next generation
		
	# Store numbers
	DisperseDeaths[gen] = sum(DisperseDeaths[gen])
	CouldNotDisperse[gen] = sum(CouldNotDisperse[gen])
	Migrants.append(len(OffDisperseIN))
	Open.append(len(tempfreegrid))
	 
	# Variables returned
	tupDoEmi = OffDisperseIN,tempfreegrid
	return tupDoEmi
	
	# End::DoEmigration()
	
# ---------------------------------------------------------------------------------------------------	
def CalculateDispersalMetrics(OffDisperseIN,xgridcopy,ygridcopy,\
Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistED,\
MDispDistED,FDispDistEDstd,MDispDistEDstd,FDispDistCD,MDispDistCD,\
FDispDistCDstd,MDispDistCDstd,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC):
	'''
	CalculateDispersalMetrics()
	This function calculates how far disperses are moving.
	'''		
	# Store the average dispersal distance offspring went
	# temp variable to store offspring dispersal distance
	FtempAvgDispDistED = []
	MtempAvgDispDistED = []
	FtempAvgDispDistCD = []
	MtempAvgDispDistCD = []
	OffDispDistCD = [] # Combined dispersal distance for storage
	Fcount = 0
	Mcount = 0
	
	# Loop through each OffDisperseIN
	for ioffspring in range(len(OffDisperseIN)):
		# Store the ED/CD distance offspring went - split up into sex
		if int(OffDisperseIN[ioffspring][0][4]) == 0:
			FtempAvgDispDistED.append(np.sqrt((xgridcopy[int(OffDisperseIN[ioffspring][0][2])]-xgridcopy[int(OffDisperseIN[ioffspring][1])])**2+(ygridcopy[int(OffDisperseIN[ioffspring][0][2])]-ygridcopy[int(OffDisperseIN[ioffspring][1])])**2))
			Fcount = Fcount + 1
			probval = Fxycdmatrix[int(OffDisperseIN[ioffspring][0][2])][int(OffDisperseIN[ioffspring][1])]
			
			# If panmictic
			if Fdispmoveno == '4' or Fdispmoveno == '6':
				cdval = 0.0
				
			# If linear
			elif Fdispmoveno == '1':
				cdval = (probval - 1.) * (-Fthreshold)
				
			# If inverse square
			elif Fdispmoveno == '2':
				
				if probval == 1.0:
					cdval = 0.0
				else:	
					cdval = np.sqrt(1. / (probval * (FScaleMax - FScaleMin) + FScaleMin))
			# If neg exponetial
			elif Fdispmoveno == '5':
				cdval = np.log((probval * (FScaleMax-FScaleMin) + FScaleMin)/float(FA)) / (-float(FB) * np.log(10))
			# If Gaussian	
			elif Fdispmoveno == '7':
				cdval = float(FB) + np.sqrt(-2*float(FC)**2 * np.log((probval*(FScaleMax-FScaleMin)+FScaleMin)/float(FA)))
			# If matrix
			elif Fdispmoveno == '8':
				cdval = (1. - probval)*(FScaleMax-FScaleMin)+FScaleMin
			elif Fdispmoveno == '9':
				cdval = probval
			
			# Write to temp	
			FtempAvgDispDistCD.append(cdval)
				
			# Then add to combined storage
			OffDispDistCD.append(cdval)
			
		# Else a male
		else:
			MtempAvgDispDistED.append(np.sqrt((xgridcopy[int(OffDisperseIN[ioffspring][0][2])]-xgridcopy[int(OffDisperseIN[ioffspring][1])])**2+(ygridcopy[int(OffDisperseIN[ioffspring][0][2])]-ygridcopy[int(OffDisperseIN[ioffspring][1])])**2))
			Mcount = Mcount + 1
			probval = Mxycdmatrix[int(OffDisperseIN[ioffspring][0][2])][int(OffDisperseIN[ioffspring][1])]
			# If panmictic
			if Mdispmoveno == '4' or Mdispmoveno == '6':
				cdval = 0.0
				MtempAvgDispDistCD.append(cdval)
			# If linear
			elif Mdispmoveno == '1':				
				cdval = (probval - 1.) * (-Mthreshold)
				MtempAvgDispDistCD.append(cdval)
			# If inverse square
			elif Mdispmoveno == '2':
				if probval == 1.0:
					cdval = 0.0
				else:	
					cdval = np.sqrt(1. / (probval * (MScaleMax - MScaleMin) + MScaleMin))
				
			# If neg exponetial
			elif Mdispmoveno == '5':
				cdval = np.log((probval * (MScaleMax-MScaleMin) + MScaleMin)/float(MA)) / (-float(MB) * np.log(10))
			elif Mdispmoveno == '7':
				cdval = float(MB) + np.sqrt(-2*float(MC)**2 * np.log((probval*(MScaleMax-MScaleMin)+MScaleMin)/float(MA)))
			elif Mdispmoveno == '8':
				cdval = (1. - probval)*(MScaleMax-MScaleMin)+MScaleMin
			elif Mdispmoveno == '9':
				cdval = probval
			MtempAvgDispDistCD.append(cdval)
			
			# Then add to combined storage
			OffDispDistCD.append(cdval)
			
	# If at least 1 Female offspring dispersed
	if Fcount > 0:		
		# And append to DispDistED
		FDispDistED.append(sum(FtempAvgDispDistED)/Fcount)
		FDispDistEDstd.append(np.std(FtempAvgDispDistED))
		# And append to DispDistCD
		FDispDistCD.append(sum(FtempAvgDispDistCD)/Fcount)
		FDispDistCDstd.append(np.std(FtempAvgDispDistCD))
	else:
		# And append to DispDistED
		FDispDistED.append(0)
		FDispDistEDstd.append(0)
		# And append to DispDistCD
		FDispDistCD.append(0)
		FDispDistCDstd.append(0)
	
	# If at least 1 Male offspring dispersed
	if Mcount > 0:		
		# And append to DispDistED
		MDispDistED.append(sum(MtempAvgDispDistED)/Mcount)
		MDispDistEDstd.append(np.std(MtempAvgDispDistED))
		# And append to DispDistCD
		MDispDistCD.append(sum(MtempAvgDispDistCD)/Mcount)
		MDispDistCDstd.append(np.std(MtempAvgDispDistCD))
	else:
		# And append to DispDistED
		MDispDistED.append(0)
		MDispDistEDstd.append(0)
		# And append to DispDistCD
		MDispDistCD.append(0)
		MDispDistCDstd.append(0)

	return OffDispDistCD
	
	# End::CalculateDispersalMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def DoDisperse(offspringno,freegrid,offspring,Fdispmoveno,Mdispmoveno,\
Fxycdmatrix,Mxycdmatrix,gen,Migrants,Open,loci,alleles,\
xgridcopy,ygridcopy,FDispDistED,MDispDistED,FDispDistCD,MDispDistCD,
logfHndl,cdevolveans,fitvals,FDispDistEDstd,MDispDistEDstd,\
FDispDistCDstd,MDispDistCDstd,subpop,subpopmigration,DisperseDeaths,CouldNotDisperse,\
gridmort,philopatry,females,subpopemigration,females_nomate,males,males_nomate,\
startSelection,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC,betas,xvars_betas,maxfit,minfit):
	'''
	DoDisperse()
	Disperse the new offspring to empty spots on grid
	Input: Units of dipsersal, movement function,
	offspring, freegrid, cdmatrix 
	Output: OffDisperseIN = [offspring,freegrid,name,[offspringgenes]] 
	'''		
		
	tupDoEmi = DoEmigration(offspring,freegrid,Migrants,Open,\
	loci,alleles,Fxycdmatrix,Mxycdmatrix,gen,offspringno,\
	cdevolveans,fitvals,subpop,subpopmigration,\
	DisperseDeaths,CouldNotDisperse,gridmort,philopatry,\
	females,subpopemigration,females_nomate,males,males_nomate,\
	startSelection,betas,xvars_betas,maxfit,minfit)
	
	OffDisperseIN = tupDoEmi[0]
	opengrids = tupDoEmi[1]
	
	# Calculate Dispersal Metrics
	OffDispDistCD = CalculateDispersalMetrics(OffDisperseIN,xgridcopy,ygridcopy,\
	Fdispmoveno,Mdispmoveno,Fxycdmatrix,Mxycdmatrix,FDispDistED,\
	MDispDistED,FDispDistEDstd,MDispDistEDstd,\
	FDispDistCD,MDispDistCD,FDispDistCDstd,MDispDistCDstd,Fthreshold,Mthreshold,FScaleMax,FScaleMin,MScaleMax,MScaleMin,FA,FB,FC,MA,MB,MC)
		
	# Return variables from this argument
	tupDoDisp = OffDisperseIN,opengrids,OffDispDistCD
	return tupDoDisp
	
	# End::DoDisperse()