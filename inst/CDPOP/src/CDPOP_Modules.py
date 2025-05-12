# -------------------------------------------------------------------------------------------------
# CDPOP_Modules.py
# Author: Erin L Landguth
# Created: June 2010
# Description: This is the function/module file for CDPOP v0.88
# --------------------------------------------------------------------------------------------------

# Import Modules with Except/Try statements

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError("Numpy required.")

# CDPOP functions
try:
	from CDPOP_PreProcess import *
except ImportError:
	raise ImportError("CDPOP PreProcess required.")
	
# Python specific functions
import os, random, copy, pdb, sys, math,itertools
from collections import Counter

# ----------------------------------------------------------
# Global symbols, if any :))
#-----------------------------------------------------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = False


# --------------------------------------------------------------------------
def PrepTextFile(textpath):
	'''
	PrepTextFile() - Prepare the input files
	'''
	
	return textpath.strip('\n').strip('\r')
	
	# End::PrepTextFile()

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
	return item,count
	
	#End::w_choice_general()
	
# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()

# ---------------------------------------------------------------------------------------------------	 
def ReadGrid(FIDnew,idnew,agenew,xgridnew,ygridnew,\
genesnew,equalsexratio,sexnew,subpopnew,infectionnew,allelst,geneswap,gen,intgenesans,hindexnew):	
	
	'''
	DoReadGrid()
	This function is reads the previous generations
	grid information at the start of the time loop.
	'''
	
	FID = FIDnew
	subpop = subpopnew
	xgrid = xgridnew
	xgridcopy = copy.deepcopy(xgridnew)
	ygrid = ygridnew
	ygridcopy = copy.deepcopy(ygridnew)
	id = idnew
	age = agenew	
	infection = infectionnew
	hindex = hindexnew
	
	# Store the number of grids
	nogrids = len(FID)
		
	# Get the number of actual filled grids and np.where indexing
	filledgrids = nogrids-sum(np.asarray(age,dtype=str)=='NA')
	
	# If gen is < geneswap, do nothing, else intialize genes here
	if (intgenesans == 'file' or intgenesans == 'random' or intgenesans == 'file_var' or intgenesans == 'random_var') and gen == geneswap:
	
		# Loop through all grids
		genes = []
		for i in range(nogrids):
									
			# And store genes information
			genes.append([])
						
			# For each loci:
			for j in range(len(allelst[0])):			
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
				for k in range(len(allelst[0][j])):
					
					# Somebody not in this spot
					if age[i] == 'NA':
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
								
						# And to genes list
						genes[i].append(tempindall)
			
			# -----------------------------
			# Update Hindex here
			# -----------------------------
			if genes[i][0] == 2:
				hindex.append(1.0)
			elif genes[i][1] == 2:
				hindex.append(0.0)
			elif genes[i][0] == 1 and genes[i][1] == 1:
				hindex.append(0.5)
			else:
				hindex.append(-9999)
	else:
		genes = genesnew
	
	# If sex is not equal ratio
	if equalsexratio == 'N' or equalsexratio == 'AtBirth':
		sex = sexnew
	
	# If equal sex ratio is Y, then split up sex into equal parts
	if equalsexratio == 'WrightFisher':
		# Error check to make sure full grids
		if sum(np.asarray(age,dtype=str)=='NA') != 0:
			print('WrightFisher is a special case assuming constant population size. Try increasing birth rate.')
			sys.exit(-1)
		sex = []
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
			
		# Then loop through each subpopulation
		for i in range(nosubpops):			
			# Then create half males and females and shuffle
			sextemp = np.append(np.zeros(subgridtotal[i]/2,"int"),np.ones(subgridtotal[i]/2,"int"))
			np.random.shuffle(sextemp)
			# Loop through these individuals and append to LIST sex
			for j in range(len(sextemp)):			
				# The add them together and shuffle and append to sex
				sex.append(str(sextemp[j]))
		
		# Delete extra stuff
		del(sextemp)
	
	#Return this functions variables
	tupReadGrid = FID,sex,id,age,xgrid,xgridcopy,ygrid,\
	ygridcopy,genes,nogrids,subpop,infection,filledgrids,hindex
	return tupReadGrid
	
	#End::DoReadGrid()

# ---------------------------------------------------------------------------------------------------	 
def GetMetrics(Population,nogrids,loci,alleles,genes,gen,Ho,\
Alleles,He,subpop,p1,p2,q1,q2,Population_age,Females_age,Males_age,age,sex,Magemort,geneswap,cdevolveans,xvars_betas,betas_selection,maxfit,minfit):
	'''
	GetMetrics()
	This function summarizes the genotypes and
	produces genetic metrics.
	Ho - Observed heterozygosity per generation
	He - Expected heterozygoisty per generation
	Alleles - Total number of unique alleles in genotype*individuals
	'''
	subpop = np.asarray(subpop)
	age = np.asarray(age)
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	if len(Magemort) > 1: # Error check
		if len(Magemort) != nosubpops:
			print('Multiple Agevars given and must match number of subpopulations specified.')
			sys.exit(-1)
	
	# Track population age numbers
	Population_age.append([]) # Add time
	Females_age.append([]) # Add time
	Males_age.append([]) # Add time
	countages = Counter(np.asarray(np.asarray(age,dtype='|U10')[np.where(np.asarray(age,dtype='|U10') != 'NA')[0]],dtype=np.int8))
	countages_F = Counter(np.asarray(np.asarray(age,dtype='|U10')[np.where(np.asarray(sex,dtype='|U10') == '0')[0]],dtype=np.int8))
	countages_M = Counter(np.asarray(np.asarray(age,dtype='|U10')[np.where(np.asarray(sex,dtype='|U10') == '1')[0]],dtype=np.int8))	
	for ipop in range(len(Magemort)): # For each subpop to track via AgeVars files, note this can be different than nosubpops
		Population_age[gen].append([])
		Females_age[gen].append([])
		Males_age[gen].append([])		
		# Get this subpop
		ipop_indexes = np.where(np.asarray(subpop,dtype='U10')==str(ipop+1))[0]		
		for iage in range(1,len(Magemort[0])+1): # For each age 1+
			Population_age[gen][ipop].append([])
			Females_age[gen][ipop].append([])
			Males_age[gen][ipop].append([])
			if len(Magemort) == 1: # For just one AgeVars file
				Population_age[gen][ipop][iage-1].append(countages[iage])
				Females_age[gen][ipop][iage-1].append(countages_F[iage])
				Males_age[gen][ipop][iage-1].append(countages_M[iage])
			else: # If more than 1 AgeVars matching population
				# Count this age in this ipop
				countage = len(np.where(np.asarray(age,dtype='|U10')[ipop_indexes]==str(iage))[0])
				Population_age[gen][ipop][iage-1].append(countage)
				countage_F = len(np.where(np.asarray(age,dtype='|U10')[ipop_indexes[np.where(np.asarray(sex,dtype='|U10')[ipop_indexes]=='0')[0]]]==str(iage))[0])
				Females_age[gen][ipop][iage-1].append(countage_F)
				countage_M = len(np.where(np.asarray(age,dtype='|U10')[ipop_indexes[np.where(np.asarray(sex,dtype='|U10')[ipop_indexes]=='1')[0]]]==str(iage))[0])
				Males_age[gen][ipop][iage-1].append(countage_M)
			
		Population_age[gen][ipop] = sum(Population_age[gen][ipop],[])
		Females_age[gen][ipop] = sum(Females_age[gen][ipop],[])
		Males_age[gen][ipop] = sum(Males_age[gen][ipop],[])
	
	# List for total, left, and right
	unique_alleles = Alleles
		
	# Only complete if greater than startGenes
	# ----------------------------------------	
	if gen >= geneswap:
		
		# ----------------------------------------------------------
		# Get summary numbers to be used in calculations or tracking		
		filledgrids = len(np.where(np.asarray(age,dtype='|U10') != 'NA')[0])# Get the number of filled grids
		
		#allele_numbers = np.asarray(range(alleles[0])*loci)# Get allele location as sequence from alleles array
		allele_numbers = []
		for iall in alleles:
			allele_numbers.append(list(range(iall)))
		allele_numbers = np.asarray(sum(allele_numbers,[]))		
		total_alleles = len(allele_numbers)	# The total number of alleles	
		
		# Get Genes array with np.nan option
		genes_array = np.asarray(genes,dtype='|U6')
		genes_array[np.where(genes_array == 'NA')[0]] = np.nan
		genes_array = np.asarray(genes_array,dtype='float')
		
		# ---------------------------------------------------
		# Get total information
		# Get allele frequency for total
		if filledgrids != 0:
			all_freq_tot = np.asarray(np.nansum(genes_array,axis=0),dtype = 'float')
			all_freq_tot = all_freq_tot/(2*filledgrids)
			#all_freq_tot = np.asarray(np.nansum(genes_array,axis=0),dtype = 'float').reshape(total_alleles)
			#all_freq_tot = all_freq_tot/(2*filledgrids)
		else:
			all_freq_tot = np.zeros(total_alleles,float)
		all_freq_sq_tot = all_freq_tot**2
		
		# Get total number of alleles
		alleles_tot = np.array(all_freq_tot>0.).sum()
		
		# Create an array to fill up with allele frequencies - only for total
		all_freq_list = np.zeros((total_alleles,2))		
		all_freq_list[:,0] = allele_numbers
		all_freq_list[:,1] = all_freq_tot
		
		# Append allele total information
		unique_alleles.append([alleles_tot])		
		
		#Calculate the number of homogenous alleles for total
		ho_count_tot = (np.array(genes_array == 2.)).sum()
				
		# Calculate the observed het for total
		if filledgrids != 0:
			ho_tot = (float(filledgrids*loci - ho_count_tot)/(loci*filledgrids))
		else:
			ho_tot = 0.0		
					
		# Append Ho information (Observed Het)
		Ho.append([ho_tot])
		
		# Calculate the homozygosity for total populations
		homozygosity_tot = sum(all_freq_sq_tot)/loci
		
		# Store He for [Total]
		if filledgrids != 0:
			he_tot = (round(1. - homozygosity_tot,2))
		else:
			he_tot = 0.0
			
		# Append He information (Expected Het)
		He.append([he_tot])
		
		# --------------------------------------------------
		# Count up individuals in each subpop and setup vars
		subgridtotal = [] # count of individuals in each subpop
		subgrids = [] # index locations
		all_freq_sub = []
		ho_count_sub = []
		ho_sub = []
		all_freq_sq_sub = []
		homozygosity_sub = []
		he_sub = []
		sumsubpopsHo = []
		sumsubpopsHe = []
		alleles_sub = []
		Population.append([]) # Tracking 
		for ipop in range(nosubpops):
			# Get location of this subpop and number in
			subgrids.append(np.asarray(np.where(np.asarray(subpop,dtype='U10') == str(ipop+1))[0][np.where(np.asarray(age,dtype='U10')[np.where(np.asarray(subpop,dtype='U10') == str(ipop+1))[0]] != 'NA')[0]],dtype='int').tolist())
			subgridtotal.append(len(np.where(np.asarray(age,dtype='U10')[np.where(np.asarray(subpop,dtype='U10') == str(ipop+1))[0]] != 'NA')[0]))
			# Add information to Population tracker
			Population[gen].append(subgridtotal[ipop])
			# for later add spots to these vars
			all_freq_sub.append([])
			ho_count_sub.append([])
			ho_sub.append([])
			all_freq_sq_sub.append([])
			homozygosity_sub.append([])
			he_sub.append([])
			alleles_sub.append([])
			
			# Allele frequency for each subpop
			if subgridtotal[ipop] != 0:
				all_freq_sub[ipop].append(np.asarray(np.sum(np.asarray(genes_array[subgrids[ipop]],dtype='float'),axis=0),dtype = 'float'))
				all_freq_sub[ipop] = all_freq_sub[ipop][0]/(2*subgridtotal[ipop])
			else:
				all_freq_sub[ipop].append(np.zeros(total_alleles,float))
				all_freq_sub[ipop] = all_freq_sub[ipop][0]
			
			# Calculate the number of homogenous alleles in each subpop
			ho_count_sub[ipop].append((np.array(np.asarray(genes_array[subgrids[ipop]],dtype='float') == 2.)).sum())
			
			# Calculate the observed het in each subpop
			if subgridtotal[ipop] != 0:
				ho_sub[ipop].append((float(subgridtotal[ipop]*loci - ho_count_sub[ipop][0])/(loci*subgridtotal[ipop])))
			else:
				ho_sub[ipop].append(0.0)
				
			# Append Ho information (Observed Het)
			Ho[gen].append(ho_sub[ipop][0])
			
			# Allele freq squared for subpops
			all_freq_sq_sub[ipop].append(all_freq_sub[ipop]**2)
			
			# Calculate the homozygosity for subpopulations
			homozygosity_sub[ipop].append(sum(all_freq_sq_sub[ipop][0])/loci)
			
			# Store He for subpopulations
			if subgridtotal[ipop] != 0:
				he_sub[ipop].append(round(1. - homozygosity_sub[ipop][0],4))
			else:
				he_sub[ipop].append(0.0)
				
			# Append He information (Expected Het)
			He[gen].append(he_sub[ipop][0])

			# Get the total number of alleles in each subpop
			alleles_sub[ipop].append(np.array(all_freq_sub[ipop]>0.).sum())
			
			# Append allele information
			unique_alleles[gen].append(alleles_sub[ipop][0])
				
		# ----------------------------------
		# Summary and fill final vars		
		Population[gen].insert(0,filledgrids) # Add Population total for tracking
				
		# Get allele frequency totals for selection section
		p1.append(all_freq_tot[0])
		p2.append(all_freq_tot[1])
		q1.append(all_freq_tot[2])
		q2.append(all_freq_tot[3])
				
		# ----------------------------------------------
		# Max/min GXE local rescaling for fitness values
		if cdevolveans.split('_')[0] == 'M':			
		
			# For Global option, at the first generation only, get all possible genotypes
			if gen == 0:
				# Calculate the total genotype space - combination replacement CR(possible alleles,2)^possible loci
				#total_genotypespace = ( math.factorial(int(cdevolveans.split('_')[3].split('A')[1]) + 2 - 1) / (math.factorial(2) * math.factorial(int(cdevolveans.split('_')[3].split('A')[1]) - 1)) )**int(cdevolveans.split('_')[2].split('L')[1])
				
				# Grab XVars as array
				xvars_betas = np.asarray(xvars_betas,dtype=float)
				if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
					multiFactor = 1
				else:
					multiFactor = 2
				
				max_linmodel2 = []
				min_linmodel2 = []
				# For each grid spot, grab the X var values
				for igrid in range(nogrids):
					grid_xvars = xvars_betas[igrid] # This grids X variables
					max_linmodel2.append([])
					min_linmodel2.append([])
					# Loop through each X var and max/min each calculation
					for ivar in range(len(grid_xvars)):
						Xvar = grid_xvars[ivar]
						betas = np.asarray(betas_selection[ivar])
						
						# if the Xvar is Zero:
						if Xvar == 0:
							# zeros out this calculation
							max_linmodel2[igrid].append(0)
							min_linmodel2[igrid].append(0)
							
						# If Xvar is positive	
						elif Xvar > 0: 
							# For max, grab max beta value (positive) at each locus and multiply by 2 or 1
							max_linmodel2[igrid].append(Xvar*np.sum(np.max(betas,1)*multiFactor))
							# For min, grab min beta value at each locus and multiply by 2 or 1
							min_linmodel2[igrid].append(Xvar*np.sum(np.min(betas,1)*multiFactor))
							
						# if Xvar is negative	
						elif Xvar < 0: 
							# For max, grab min beta value at each locus and multiply by 2 or 1
							max_linmodel2[igrid].append(Xvar*np.sum(np.min(betas,1)*multiFactor))
							# For min, grab max beta value (positive) at each locus and multiply by 2 or 1 
							min_linmodel2[igrid].append(Xvar*np.sum(np.max(betas,1)*multiFactor))	
						
					max_linmodel2[igrid] = sum(max_linmodel2[igrid])
					min_linmodel2[igrid] = sum(min_linmodel2[igrid])
				
				maxfit.append(max(max_linmodel2))
				minfit.append(min(min_linmodel2))		
				'''				
				# Get possible allele combinations with replacement, but not duplicates
				range_alleles = range(int(cdevolveans.split('_')[3].split('A')[1])) 
				possible_alleles_index = list(itertools.combinations_with_replacement(range_alleles, 2))
						
				# The translate the above to 2/1/0
				possible_alleles_list = []
				for i in xrange(len(possible_alleles_index)):
					possible_alleles_list.append(np.zeros(int(cdevolveans.split('_')[3].split('A')[1])).tolist()) 
					for thisall in possible_alleles_index[i]:
						possible_alleles_list[i][thisall] = possible_alleles_list[i][thisall] + 1
				
				# Get cartesian product, repeat number of loci times 
				possible_genes = list(itertools.product(possible_alleles_list,repeat=int(cdevolveans.split('_')[2].split('L')[1])))
				
				# Calculate the linmodel for the possible alleles combinations - max/min values
				# -----------------------------------------------
				#pdb.set_trace()
				max_linmodel = []
				min_linmodel = []
				xvars_betas = np.asarray(xvars_betas,dtype=float)
				# Loop through each individual spot
				for iind in xrange(len(possible_genes)):
					selgenes = sum(possible_genes[iind],[])
					max_linmodel.append([]) # Add spot in array
					min_linmodel.append([]) # Add spot in array
					# Loop through the environmental vars
					for ixvar in xrange(int(cdevolveans.split('_')[1].split('X')[1])):  
						# Get max and min xvar in this vector
						min_xvar = min(xvars_betas[:,ixvar])
						max_xvar = max(xvars_betas[:,ixvar])
						# loop through the alleles
						betas = sum(betas_selection[ixvar],[])
						for iall in xrange(len(selgenes)):
							thisbeta = betas[iall]
							thisallele = selgenes[iall]
							if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
								if thisallele == 2:
									thisallele = 1 # Change second copy to 1
							# allele X environment X beta effect
							val1 = min_xvar * thisbeta * thisallele # because of signs on betas and xvars, this valuues might not be max and mins, so check
							val2 = max_xvar * thisbeta * thisallele
							if val1 <= val2: # minimum value
								min_linmodel[iind].append(val1)
								max_linmodel[iind].append(val2)
							else:
								min_linmodel[iind].append(val2)
								max_linmodel[iind].append(val1)
					
					# Get Epistasis term if applicable
					if epistasis != 'N':
						print('Epistasis still in beta.')
						sys.exit(-1)
						epistasis_beta = float(epistasis.split('_')[0])
						epistasis_Alocispot = int(epistasis.split('_')[1].split('A')[1][0])-1 # count from zero subtract 1
						epistasis_Aallspot = int(epistasis.split('_')[1].split('A')[1][1])-1 # count from zero subtract 1
						epistasis_Agene = selgenes[epistasis_Alocispot*alleles[0]+epistasis_Aallspot]
						epistasis_Blocispot = int(epistasis.split('_')[2].split('B')[1][0])-1 # count from zero subtract 1
						epistasis_Ballspot = int(epistasis.split('_')[2].split('B')[1][1])-1 # count from zero subtract 1
						epistasis_Bgene = selgenes[epistasis_Blocispot*alleles[0]+epistasis_Ballspot]			
						
					# Add the beta not and sum
					min_linmodel[iind].append(betas_selection[-1])
					min_linmodel[iind] = sum(min_linmodel[iind])
					max_linmodel[iind].append(betas_selection[-1])
					max_linmodel[iind] = sum(max_linmodel[iind])
				
				# Get the max/min GXE calculation
				maxfit.append(max(max_linmodel))
				minfit.append(min(min_linmodel))
				pdb.set_trace()
				del(possible_genes)				
			
			# Code section for 'local option' to calculate linmodel at each gen
			# Get linmodel value for all individuals - 
			# --------------------------------------
			# Get selection-driven region		
			indexto = int(cdevolveans.split('_')[2].split('L')[1]) * int(cdevolveans.split('_')[3].split('A')[1])
			selgenes = genes_array[:,0:indexto]	
			linmodel = []
			# Loop through each individual spot
			for iind in xrange(len(selgenes)): 
				# If there is an NA here
				if selgenes[iind][0] == 'NA':
					pdb.set_trace()
				linmodel.append([]) # Add spot in array
				# Loop through the environmental vars
				for ixvar in xrange(int(cdevolveans.split('_')[1].split('X')[1])):  
					xvar = float(xvars_betas[iind][ixvar])
					# Loop through the alleles
					betas = sum(betas_selection[ixvar],[])
					for iall in xrange(len(selgenes[iind])):
						thisbeta = betas[iall]
						thisallele = selgenes[iind][iall]
						if cdevolveans.split('_')[4] == 'ModelY': # Code 1,0
							if thisallele == 2:
								thisallele = 1 # Change second copy to 1
						# allele X environment X beta effect
						linmodel[iind].append(xvar * thisbeta * thisallele)			
				# Add the beta not and sum
				linmodel[iind].append(betas_selection[-1])
				linmodel[iind] = sum(linmodel[iind])
			
			# Get the max/min GXE calculation
			maxfit.append(max(linmodel))
			minfit.append(min(linmodel))		
			'''

	# Return and store numbers if less than geneswap time
	else:
		# --------------------------------------------------
		# Count up individuals in each subpop and setup vars
		subpop = np.asarray(subpop)
		age = np.asarray(age)
		subgridtotal = [] # count of individuals in each subpop
		Population.append([]) # Tracking
		p1.append(np.nan)
		p2.append(np.nan)
		q1.append(np.nan)
		q2.append(np.nan)
		unique_alleles.append([np.nan])
		He.append([np.nan])
		Ho.append([np.nan])
		for ipop in range(nosubpops):
			subgridtotal.append(len(np.where(np.asarray(age,dtype='U10')[np.where(np.asarray(subpop,dtype='U10') == str(ipop+1))[0]] != 'NA')[0]))
			# Add information to Population tracker
			Population[gen].append(subgridtotal[ipop])
			unique_alleles[gen].append(np.nan)
			He[gen].append(np.nan)
			Ho[gen].append(np.nan)
		# Add one more for total
		unique_alleles[gen].append(np.nan)
		He[gen].append(np.nan)
		Ho[gen].append(np.nan)		
		
		# And then get the number of filled grids
		filledgrids = len(np.where(np.asarray(age,dtype='U10') != 'NA')[0])
		# Add Population total
		Population[gen].insert(0,filledgrids)
			
	subpop = subpop.tolist()
	age = age.tolist()
	# Return variables from this function
	tupGetMetrics = filledgrids,subgridtotal
	return tupGetMetrics
	
	#End::GetMetrics()
	
# ---------------------------------------------------------------------------------------------------	 
def InheritGenes(gen,AllelesMutated,offspringno,\
offspring,genes,loci,muterate,mtdna,mutationans,geneswap,epireset,Track_EpigeneReset1,Track_EpigeneReset2,startEpigene,epigeneans,cdevolveans,noalleles,hindex):
	'''
	InheritGenes()
	Pass along gentic information to survived offspring from parents
	Input: offspring, genes 
	Output: [femaleid,maleid,cdmatidofmother,cdmatidoffather,sex,infection,TWINID,[
	genetic information], hindex]		
	'''		
	
	Track_EpigeneReset1.append([]) # Tracking for resets
	Track_EpigeneReset2.append([]) # Tracking for resets
	# Get epiloci region
	if epigeneans != 'N':
		if cdevolveans.split('_')[0] == 'M':
			# Then the first l loci are for selection, and next epigenetics
			selloci = int(cdevolveans.split('_')[2].split('L')[1])
			epiloci = int(epigeneans.split('_')[1].split('L')[1])
		else:
			selloci = 0
			epiloci = int(epigeneans.split('_')[1].split('L')[1])
		epiloci_index = list(range(selloci*int(cdevolveans.split('_')[3].split('A')[1]),selloci*int(cdevolveans.split('_')[3].split('A')[1])+epiloci*int(epigeneans.split('_')[2].split('A')[1]))) # Index location for epi region, epi region only considers 2 alleles but generic form to index to here.
	else:
		epiloci_index = [-9999]
	#pdb.set_trace()
	# Check for generation to start swapping genes
	if gen >= geneswap:
		
		# Storage for tracking how many alleles mutated
		noallelesmutated = []			
		# If there are offspring
		if int(offspringno) != int(0):
			
			# Sort by twins. 
			offspring.sort(key=lambda x: x[6]) 
			tempoff = np.asarray(copy.deepcopy(offspring))
			countTwins = Counter(tempoff[:,6])
	
			isTwin = False # For checking for twins
			twingenes = [] # Initialize twin genes
			# Begin loop through offspring
			for i in range(offspringno):				
				
				# If twins genes were copied already from previous offspring, skip this section
				if len(twingenes) == 0: 
					# Temp storage for i's mother's  and fathers genes
					mothergenes=genes[int(offspring[i][0])]
					fathergenes=genes[int(offspring[i][1])]					
					fathergenes = np.asarray(fathergenes,dtype=int)
					mothergenes = np.asarray(mothergenes,dtype=int)
					# Temp genes storage for offspring
					tempgenes = np.zeros(len(fathergenes),dtype =int)
					# Allele indices
					alleles = np.asarray(list(range(len(mothergenes))))
					
					# Loop through loci
					for iloci in range(loci): 
						# Allele indices to sample from - index into tempgenes
						#possiblealleles = alleles[(iloci*len(mothergenes)/loci):(iloci*len(mothergenes)/loci+len(mothergenes)/loci)]
						possiblealleles = alleles[sum(noalleles[0:iloci]):sum(noalleles[0:iloci+1])]
						#pdb.set_trace() # if statement not right with new epiloci check	
						# If this is not the epiregion, assume diploid, randomly grab from parents alleles
						# ----------------------------------------------------------------
						if len(np.where(np.asarray(epiloci_index) == iloci*2)[0]) == 0:
						#if epiloci_index != range(iloci*2,iloci*2+epiloci*2):
						# Is this loci part of the epiloci region and index? *2 for the new flattened epiloci index 2 alleles. 
							# Father and mother locations							
							F2 = np.where(fathergenes[possiblealleles] == 2)[0] # location of 2s
							F1 = np.where(fathergenes[possiblealleles] == 1)[0]
							M2 = np.where(mothergenes[possiblealleles] == 2)[0]
							M1 = np.where(mothergenes[possiblealleles] == 1)[0]
							Falls = np.concatenate((F2,F2,F1),axis=0) # 2 copies of 2s
							Malls = np.concatenate((M2,M2,M1),axis=0) # 2 copies of 2s		
							# Sample allele from each parent
							FsampleAlleles = np.random.choice(Falls,1).tolist()
							MsampleAlleles = np.random.choice(Malls,1).tolist()
							# Fill in alleles corresponding to sampled spots
							tempgenes[possiblealleles[FsampleAlleles[0]]] = tempgenes[possiblealleles[FsampleAlleles[0]]] + 1
							tempgenes[possiblealleles[MsampleAlleles[0]]] = tempgenes[possiblealleles[MsampleAlleles[0]]] + 1
						
						# Epiregion, not necessarily diploid, different checks, and also check for resets
						# --------------------------------------------------------------
						else:
							pdb.set_trace() # Check section for flattened gene array
														
							# Get reset numbers for each allele
							Reset1 = float(epireset[offspring[i][0]][np.where(np.asarray(epiloci_index)==iloci*2)[0][0]].split(';')[0]) # first allele reset
							Reset2 = float(epireset[offspring[i][0]][np.where(np.asarray(epiloci_index)==iloci*2)[0][0]].split(';')[1]) # second allele reset
							
							# Randomly grab one of the alleles from each parent
							Fsample_allindex = np.random.choice(possiblealleles,1).tolist() # father's index location
							Msample_allindex = np.random.choice(possiblealleles,1).tolist() # mother's index location
							
							# Offspring inherits these epialleles or not from parents
							Fsample_allval = fathergenes[Fsample_allindex[0]] # father's allele value
							Msample_allval = mothergenes[Msample_allindex[0]] # mother's allele value
							# Check homo cases, offspring only gets one copy from each parent
							if Fsample_allval == 2:
								Fsample_allval = 1
							if Msample_allval == 2:
								Msample_allval = 1
							# Fill in offspring genes, indexing to allele location
							tempgenes[Fsample_allindex[0]] = tempgenes[Fsample_allindex[0]] + Fsample_allval
							tempgenes[Msample_allindex[0]] = tempgenes[Msample_allindex[0]] + Msample_allval
							
							# Check for resets
							# ----------------
							# Independent alleles
							if epigeneans.split('_')[4] == 'Ind':
								# First allele
								# ------------
								if tempgenes[possiblealleles[0]] == 2: # two copies, two reset checks
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # This allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # This allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
								elif tempgenes[possiblealleles[0]] == 1: # one copy, one reset check
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # This allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
								else: # 0 copies, pass
									pass
								# Second allele
								# -------------
								if tempgenes[possiblealleles[1]] == 2: # two copies, two reset checks
									rand_reset = np.random.uniform()							
									if rand_reset < Reset2: # This allele resets
										tempgenes[possiblealleles[1]] = tempgenes[possiblealleles[1]] - 1
										Track_EpigeneReset2[gen].append(1)
									rand_reset = np.random.uniform()							
									if rand_reset < Reset2: # This allele resets
										tempgenes[possiblealleles[1]] = tempgenes[possiblealleles[1]] - 1
										Track_EpigeneReset2[gen].append(1)
								elif tempgenes[possiblealleles[1]] == 1: # one copy, one reset check
									rand_reset = np.random.uniform()							
									if rand_reset < Reset2: # This allele resets
										tempgenes[possiblealleles[1]] = tempgenes[possiblealleles[1]] - 1
										Track_EpigeneReset2[gen].append(1)
								else: # 0 copies, pass
									pass
							
							# Methylation, dependent case
							# ---------------------------
							else:
								# First allele
								# ------------
								if tempgenes[possiblealleles[0]] == 2: # two copies, two reset checks
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # This allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # This allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
									# Second allele
									# -------------
									if tempgenes[possiblealleles[1]] != 0:
										# Should not have a copy of this allele
										print('Error in epigenetic locus.')
										sys.exit(-1)
								elif tempgenes[possiblealleles[0]] == 1: # one copy, one reset check
									rand_reset = np.random.uniform()							
									if rand_reset < Reset1: # First allele resets
										tempgenes[possiblealleles[0]] = tempgenes[possiblealleles[0]] - 1
										Track_EpigeneReset1[gen].append(1)
										# Second allele
										# ------------
										if tempgenes[possiblealleles[1]] == 2: 
											# Should not have a 2 copies of this second allele
											print('Error in epigenetic locus.')
											sys.exit(-1)
										elif tempgenes[possiblealleles[1]] == 1:
											# If first allele reset, then second one automatically resets
											tempgenes[possiblealleles[1]] = 0
											Track_EpigeneReset2[gen].append(1)
										else: # Does not have a copy of this second allele, pass
											pass
									else: # First allele did not reset
										# Check Second allele, 
										# --------------------
										if tempgenes[possiblealleles[1]] == 2: 
											# Should not have a 2 copies of this allele
											print('Error in epigenetic locus.')
											sys.exit(-1)
										elif tempgenes[possiblealleles[1]] == 1:
											# Possibility that methylation not inherited
											rand_reset = np.random.uniform()							
											if rand_reset < Reset2: # Second allele resets
												tempgenes[possiblealleles[1]] = 0
												Track_EpigeneReset2[gen].append(1)
										else: # Does not have a copy of this second allele, pass
											pass
								else:	# 0 copies of first allele
									# Check second allele
									# -------------------
									if tempgenes[possiblealleles[1]] == 2:
										# Assume that second allele not passed on either
										tempgenes[possiblealleles[1]] = 0
										Track_EpigeneReset2[gen].append(2)
									elif tempgenes[possiblealleles[1]] == 1:
										# Assume that second allele not passed on either
										tempgenes[possiblealleles[1]] = 0
										Track_EpigeneReset2[gen].append(1)
									else: #0 copies, pass
										pass		
												
					# mtDNA is turned on
					# ------------------
					if mtdna == 'Y':
						# Force last locus to be mothergenes - possible alleles are from the last loop above
						tempgenes[possiblealleles] = mothergenes[possiblealleles]				
						
					# Check for twin copy
					# -------------------
					if len(offspring[i][6].split('T')) > 1:
						# Check if there are indeed more than 1 of these 
						if countTwins[offspring[i][6]] > 1:
							twingenes = copy.deepcopy(tempgenes)
										
				# Then flip the twin check back for next offsprings
				else: # This is a twin
					tempgenes = copy.deepcopy(twingenes)
					isTwin = True
					# Check
					if offspring[i][6] != offspring[i-1][6]:
						print('Twinning algorithm is not working. Email Erin.')
						sys.exit(-1)
														
				# Then check for mutations at each allele
				# ---------------------------------------
				if muterate != 0.0:
					for iloci in range(loci): # Loop through loci
						mutationrandnos = np.random.uniform(size=2) # Get a random number for checking
						# Allele indices to sample from - index into tempgenes
						#possiblealleles = alleles[(iloci*len(mothergenes)/loci):(iloci*len(mothergenes)/loci+len(mothergenes)/loci)]
						possiblealleles = alleles[sum(noalleles[0:iloci]):sum(noalleles[0:iloci+1])] # This accounts for variable alleles per loci.
						# Get the current location of alleles - index into tempgenes 
						thisloci = possiblealleles[np.where(tempgenes[possiblealleles] != 0)[0]]
						# Check case for homo
						if len(thisloci) == 1:
							# Copy the spot
							thisloci = np.concatenate((thisloci,thisloci),axis=0)
						
						# Loop through alleles
						for iall in range(2): 			
							
							# Check if random number is less than muterate
							if mutationrandnos[iall] < muterate:
								
								# First remove this allele from tempgenes
								tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] - 1
																
								# If random kth allele model
								if mutationans == 'random':
									# Randomly choose another allele, but not what allele it was									
									movealleleTO = np.random.choice(possiblealleles[np.where(thisloci[iall] != possiblealleles)[0]],1).tolist()[0]
									# Index into tempgenes and add 1
									tempgenes[movealleleTO] = tempgenes[movealleleTO] + 1
																	
									# Count a mutation
									noallelesmutated.append(1)
									
								# If just forward mutation
								elif mutationans == 'forward':
									# Move allele forward unless it is the last one
									if thisloci[iall] != possiblealleles[-1]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)
									else:	
										# The allele did not mutate because it was on the ends
										tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] + 1

								# If just forward mutation
								elif mutationans == 'backward':
									# Move allele backward unless it is the first one
									if thisloci[iall] != possiblealleles[0]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)
									else:	
										# The allele did not mutate because it was on the ends
										tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] + 1
										
								# If forward and backward mutation
								elif mutationans == 'forwardbackward':
									# Then random forward or backward step
									randstep = np.random.uniform()
									# To go left, but it can't be the first allele
									if randstep < 0.5 and thisloci[iall] != possiblealleles[0]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
										# Count a mutation
										noallelesmutated.append(1)
									# To go right, but it can't be the last allele
									elif randstep >= 0.5 and thisloci[iall] != possiblealleles[-1]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1
										# Count a mutation
										noallelesmutated.append(1)
									else:	
										# The allele did not mutate because it was on the ends
										tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] + 1
								
								# If forward mutation in A and backward mutation for b (A -> a, b -> B)
								elif mutationans == 'forwardAbackwardBrandomN':
									print('Currently not operating. Email Erin.')
									sys.exit(-1)
									if iloci == 0 and thisloci[iall] == possiblealleles[0]:
										tempgenes[thisloci[iall]+1] = tempgenes[thisloci[iall]+1] + 1									
										# Count a mutation
										noallelesmutated.append(1)
									elif iloci == 1 and thisloci[iall] == possiblealleles[1]:
										tempgenes[thisloci[iall]-1] = tempgenes[thisloci[iall]-1] + 1
										# Count a mutation
										noallelesmutated.append(1)
									elif iloci != 0 and iloci != 1:
										# Randomly choose another allele								
										movealleleTO = np.random.choice(possiblealleles[np.where(thisloci[iall] != possiblealleles)[0]],1).tolist()[0]
										# Index into tempgenes and add 1
										tempgenes[movealleleTO] = tempgenes[movealleleTO] + 1
																		
										# Count a mutation
										noallelesmutated.append(1)
									else:	
										# The allele did not mutate because it was on the ends
										tempgenes[thisloci[iall]] = tempgenes[thisloci[iall]] + 1
								
								# No other mutation models matched
								else:
									print('The mutation model does not exist.')
									sys.exit(-1)	
				
				# Add to offspring list
				# ---------------------
				# Add spot in offspring array for individual i's genes
				offspring[i].append(tempgenes.tolist())
				
				# Add Hindex to offspring list
				# ----------------------------
				M_hindex = float(hindex[int(offspring[i][0])])
				F_hindex = float(hindex[int(offspring[i][1])])
				off_hindex = M_hindex/2. + F_hindex/2.
				offspring[i].append(off_hindex)				
				
				# Copy the saved twingenes for the next offspring
				if isTwin:
					twingenes = [] # Erase the copied twin genes
					isTwin = False # Turn off 
					
			# Now store the total number of alleles that mutated
			AllelesMutated.append(sum(noallelesmutated))
			
			# Delete temp variables to free up space
			del(mothergenes)
			del(fathergenes)
			#del(tempgenes)
			
		# If there are no offspring
		elif int(offspringno) == int(0):
		
			# Store the total number of alleles that mutated
			AllelesMutated.append(0.0)				
	
	# If not generation to swap genes
	else:
		
		# Store the total number of alleles that mutated
		AllelesMutated.append(0.0)	
		
		# Add Hindex to offspring list
		
	
	# Tracking numbers, summary	
	Track_EpigeneReset1[gen] = sum(Track_EpigeneReset1[gen])
	Track_EpigeneReset2[gen] = sum(Track_EpigeneReset2[gen])
		
	# Return variables from this argument
	return offspring
	
	# End::InheritGenes()

# ---------------------------------------------------------------------------------------------------
def ConstantMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,FID,Magemort,Fagemort,infection,geneswap,mature,subpop,hindex):
	
	'''
	Constant mortality applied to each age class, track total number of deaths for each 
	age and sex class.
	'''	
	
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Split up for sex
	females = np.where(np.asarray(sex,dtype='|U4')=='0')[0]
	males = np.where(np.asarray(sex,dtype='|U4')=='1')[0]
	# Split for subpops
	Fpops = np.asarray(subpop)[females]
	Mpops = np.asarray(subpop)[males]
		
	Fagedeaths = []
	Magedeaths = []
	extra_Fagedeaths = [] # For older than # in age class
	extra_Magedeaths = [] # For older than # in age class
	Mdeleteoldindex = []	
	Fdeleteoldindex = []
	# Loop through subpops, but only if more than 1 AgeVars file is given
	if len(Magemort) > 1:			
		for ipop in range(len(Magemort)):
			Fagedeaths.append([])
			Magedeaths.append([])
			extra_Fagedeaths.append([])
			extra_Magedeaths.append([])
			Mdeleteoldindex.append([])	
			Fdeleteoldindex.append([])
			# Get all females/males in this pop - get original index into females/males
			females_inthispop = females[np.where(Fpops == str(ipop+1))[0]]
			males_inthispop = males[np.where(Mpops == str(ipop+1))[0]]
			
			# First females: If there are females in this subpop
			if len(females_inthispop) != 0:
				# Get unique ages
				Fages = np.asarray(age)[females_inthispop]				
				Funiqueages = Counter(np.asarray(np.asarray(Fages,dtype='|U10')[np.where(np.asarray(Fages,dtype='|U10') != 'NA')[0]],dtype=np.int8))
								
				# Loop through unique ages, counting deaths, age 1+ incorporated here
				for iage in range(1,(len(Fagemort[ipop])+1)):
					Fagedeaths[ipop].append(int(round(Fagemort[ipop][iage-1]*Funiqueages[iage])))
					# Then take sample to delete from 
					Fdeleteoldindex[ipop].append(np.random.choice(females_inthispop[np.where(np.asarray(Fages,dtype = 'str') == str(iage))[0]],Fagedeaths[ipop][iage-1],replace=False).tolist())
								
				# Check for individuals older than age class number
				if max(Funiqueages) > len(Fagemort[ipop]):
					print('Warning: Female age class exceeds specified class in Agevars.csv file. Recommend 100% mortality for last age class. Grouping these age classes.')
					for j in range(len(Fagemort[ipop])+1,(max(Funiqueages)+1)):
						extra_Fagedeaths[ipop].append(int(round(Fagemort[ipop][-1]*Funiqueages[j])))
					# If there are extra deaths then sample to delete from
					if len(extra_Fagedeaths[ipop]) != 0:
						count = len(Fagemort[ipop])+1
						for k in range(len(extra_Fagedeaths[ipop])):
							Fdeleteoldindex[ipop].append(np.random.choice(females_inthispop[np.where(np.asarray(Fages,dtype='|U10') == str(count))[0]],extra_Fagedeaths[ipop][k],replace=False).tolist())
							count = count + 1
					
			# Then males: If there are males in this subpop
			if len(males_inthispop) != 0:
				# Get unique ages
				Mages = np.asarray(age)[males_inthispop]
				Muniqueages = Counter(np.asarray(np.asarray(Mages,dtype='|U10')[np.where(np.asarray(Mages,dtype='|U10') != 'NA')[0]],dtype=np.int8))
				
				# Loop through unique ages, counting deaths, age 1+ incorporated here
				for iage in range(1,(len(Magemort[ipop])+1)):
					Magedeaths[ipop].append(int(round(Magemort[ipop][iage-1]*Muniqueages[iage])))
					# Then take sample to delete from 
					Mdeleteoldindex[ipop].append(np.random.choice(males_inthispop[np.where(np.asarray(Mages,dtype = 'str') == str(iage))[0]],Magedeaths[ipop][iage-1],replace=False).tolist())
				
				# Check for individuals older than age class number
				if max(Muniqueages) > len(Magemort[ipop]):
					print('Warning: Male age class exceeds specified class in Agevars.csv file. Recommend 100% mortality for last age class. Grouping these age classes.')
					for j in range(len(Magemort[ipop])+1,(max(Muniqueages)+1)):
						extra_Magedeaths[ipop].append(int(round(Magemort[ipop][-1]*Muniqueages[j])))
					# If there are extra deaths then sample to delete from
					if len(extra_Magedeaths[ipop]) != 0:
						count = len(Magemort[ipop])+1
						for k in range(len(extra_Magedeaths[ipop])):
							Mdeleteoldindex[ipop].append(np.random.choice(males_inthispop[np.where(np.asarray(Mages,dtype='|U10') == str(count))[0]],extra_Magedeaths[ipop][k],replace=False).tolist())
							count = count + 1
	# If there was just one AgeVars file given
	else:
		ipop = 0
		Fagedeaths.append([])
		Magedeaths.append([])
		extra_Fagedeaths.append([])
		extra_Magedeaths.append([])
		Mdeleteoldindex.append([])	
		Fdeleteoldindex.append([])
				
		# First females: If there are females in this subpop
		if len(females) != 0:
			# Get unique ages
			Fages = np.asarray(age)[females]
			Funiqueages = Counter(np.asarray(np.asarray(Fages,dtype='|U10')[np.where(np.asarray(Fages,dtype='|U10') != 'NA')[0]],dtype=np.int8))
			# Loop through unique ages, counting deaths, age 1+ incorporated here
			for iage in range(1,(len(Fagemort[ipop])+1)):
				Fagedeaths[ipop].append(int(round(Fagemort[ipop][iage-1]*Funiqueages[iage])))
				# Then take sample to delete from 
				Fdeleteoldindex[ipop].append(np.random.choice(females[np.where(np.asarray(Fages,dtype = 'str') == str(iage))[0]],Fagedeaths[ipop][iage-1],replace=False).tolist())
						
			# Check for individuals older than age class number
			if max(Funiqueages) > len(Fagemort[ipop]):
				print('Warning: Female age class exceeds specified class in Agevars.csv file. Recommend 100% mortality for last age class. Grouping these age classes.')
				for j in range(len(Fagemort[ipop])+1,(max(Funiqueages)+1)):
					extra_Fagedeaths[ipop].append(int(round(Fagemort[ipop][-1]*Funiqueages[j])))
				# If there are extra deaths then sample to delete from
				if len(extra_Fagedeaths[ipop]) != 0:
					count = len(Fagemort[ipop])+1
					for k in range(len(extra_Fagedeaths[ipop])):
						Fdeleteoldindex[ipop].append(np.random.choice(females[np.where(np.asarray(Fages,dtype='|U10') == str(count))[0]],extra_Fagedeaths[ipop][k],replace=False).tolist())
						count = count + 1
				
		# Then males: If there are males in this subpop
		if len(males) != 0:
			# Get unique ages
			Mages = np.asarray(age)[males]
			Muniqueages = Counter(np.asarray(np.asarray(Mages,dtype='|U10')[np.where(np.asarray(Mages,dtype='|U10') != 'NA')[0]],dtype=np.int8))
			# Loop through unique ages, counting deaths, age 1+ incorporated here
			for iage in range(1,(len(Magemort[ipop])+1)):
				Magedeaths[ipop].append(int(round(Magemort[ipop][iage-1]*Muniqueages[iage])))
				# Then take sample to delete from 
				Mdeleteoldindex[ipop].append(np.random.choice(males[np.where(np.asarray(Mages,dtype = 'str') == str(iage))[0]],Magedeaths[ipop][iage-1],replace=False).tolist())
			
			# Check for individuals older than age class number
			if max(Muniqueages) > len(Magemort[ipop]):
				print('Warning: Male age class exceeds specified class in Agevars.csv file. Recommend 100% mortality for last age class. Grouping these age classes.')
				for j in range(len(Magemort[ipop])+1,(max(Muniqueages)+1)):
					extra_Magedeaths[ipop].append(int(round(Magemort[ipop][-1]*Muniqueages[j])))
				# If there are extra deaths then sample to delete from
				if len(extra_Magedeaths[ipop]) != 0:
					count = len(Magemort[ipop])+1
					for k in range(len(extra_Magedeaths[ipop])):
						Mdeleteoldindex[ipop].append(np.random.choice(males[np.where(np.asarray(Mages,dtype='|U10') == str(count))[0]],extra_Magedeaths[ipop][k],replace=False).tolist())
						count = count + 1
	
	# Flatten and turn into array and get original index locations and add indices together
	if sum(sum(Mdeleteoldindex,[]),[]) != []:
		Mdeleteoldindex = np.asarray(sum(sum(Mdeleteoldindex,[]),[]))
		deleteallindex = np.append(openindex,Mdeleteoldindex)
	else:
		Mdeleteoldindex = np.asarray([])
		deleteallindex = openindex
	if sum(sum(Fdeleteoldindex,[]),[]) != []:
		Fdeleteoldindex = np.asarray(sum(sum(Fdeleteoldindex,[]),[]))
		deleteallindex = np.append(deleteallindex,Fdeleteoldindex)
	else:
		Fdeleteoldindex = np.asarray([])
	deleteallindex = np.asarray(deleteallindex,dtype = 'int')
	
	# Store freegrid locations
	freegrid = np.asarray(FID)[deleteallindex]
	freegrid = list(freegrid)	
	# Delete all of the old generations information.
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	hindex = np.delete(np.asarray(hindex),deleteallindex)
	
	# Keep genes around if within burnin faze
	if gen >= geneswap:
		genes = np.delete(genes,deleteallindex,axis=0)
		#genes = np.asarray(genes,dtype=int)
		genes = genes.tolist()
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)
	
	# Store total number of Deaths information for output	
	Track_MDeaths.append(Magedeaths)
	Track_FDeaths.append(Fagedeaths)
	# Check that agedeaths equals age distribution file
	# Commenting this out for now...
	#if len(extra_Magedeaths) != 0:
	#	Track_MDeaths[gen][-1] = Track_MDeaths[gen][-1] + sum(extra_Magedeaths)
	#if len(extra_Fagedeaths) != 0:
	#	Track_FDeaths[gen][-1] = Track_FDeaths[gen][-1] + sum(extra_Fagedeaths)
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,mature,hindex
	return tupMort
	# End::DoConstantMortality
	
# ----------------------------------------------------------------------------------------------	 
def DDMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,FID,infection,geneswap,K_env,popmodel,Magemort,Fagemort,mature,subpop,hindex):
	'''
	DensityDependentMortality()
	Density dependent survival applied to each population.		
	'''
	
	# Multiple Age Vars files check
	if len(Magemort) > 1:
		print('Multiple Agevars files given, and DDmortality option not currently operating.')
		sys.exit(-1)
	
	# Get total number of deaths for each age class, be careful of NAs and strings
	uniqueages = Counter(np.asarray(np.asarray(age,dtype='|U10')[np.where(np.asarray(age,dtype='|U10') != 'NA')[0]],dtype=np.int8))
	Nt = len(np.where(np.asarray(age,dtype='|U10') != 'NA')[0])
	agedeaths = []	
	extra_agedeaths = []
	
	if Magemort[0] != Fagemort[0]:
		print('Warning: Logistic growth is specified and the average age specific mortality of males and females will be used.')
	agemort = (np.asarray(Magemort[0]) + np.asarray(Fagemort[0])) / 2.
	
	for i in range(1,(len(agemort)+1)): # Only for ages 1 +		
		# Current Ni,t
		Ni = uniqueages[i]		
		# If last age class, no survivors
		if i == len(agemort):			
			mortreturn = Ni
		else:
			# Next age class
			Niplus1 = uniqueages[i+1]
			
			# Rickers option
			if popmodel == 'Rickers' or popmodel == 'rickers':
				print('Not operational.')
				sys.exit(-1)
			# Richards
			elif popmodel == 'Richards' or popmodel == 'richards':
				print('Not operational.')
				sys.exit(-1)
			# Logistic
			elif popmodel == 'logistic' or popmodel == 'Logistic':
				# Ni+1,t * (1 - Nt/K) * (si-1 * Ni,t - Ni+1,t)
				Ntplus1 = Niplus1 + (1. - Nt/float(K_env))*((1.-agemort[i-1]) * Ni - float(Niplus1))
				mortreturn = Ni - int(Ntplus1)
				if mortreturn < 0:
					mortreturn = 0
			else:
				print('Enter rickers, richards, or logistic for density dependent option')
				sys.exit(-1)
		
		agedeaths.append(mortreturn)
	
	# Switch for ages over age classes, apply last mortality to age classes
	if max(uniqueages) > len(agemort):
		print('Warning: age classes exceeding specified class in Agevars.csv file, apply 0 survival.')
		for j in range(len(agemort)+1,(max(uniqueages)+1)):
			extra_agedeaths.append(uniqueages[j])
			
	# Grab locations that are open
	openindex = np.where(np.asarray(sex) == 'NA')[0]
	
	# Grab locations that are not open
	filledindex = np.where(np.asarray(sex) != 'NA')[0]
	
	# Then take a sample from the possible age class indices to delete from
	deleteoldindex = []
	for i in range(1,(len(agedeaths)+1)): # index from 1 to ...
		# NA switch
		if len(openindex) == 0:
			deleteoldindex.append(np.random.choice(np.where(np.asarray(age) == i)[0],int(agedeaths[i-1]),replace=False).tolist()) # index into age deaths - 1
		else:
			deleteoldindex.append(np.random.choice(np.where(np.asarray(age) == str(i))[0],int(agedeaths[i-1]),replace=False).tolist())	
	# In case there are extra age deaths
	if len(extra_agedeaths) != 0:
		count = len(agemort)+1
		for j in range(len(extra_agedeaths)):
			deleteoldindex.append(np.random.choice(np.where(np.asarray(age) == count)[0],int(extra_agedeaths[j]),replace=False).tolist())
			count = count + 1
	
	# Flatten and turn into array
	deleteoldindex = np.asarray([item for sublist in deleteoldindex for item in sublist],dtype = 'int')
		
	# Then add indices together
	deleteallindex = np.append(openindex,deleteoldindex)
	deleteallindex = np.asarray(deleteallindex,dtype='int')
	
	# Store freegrid locations
	freegrid = np.asarray(FID)[deleteallindex]
	freegrid = list(freegrid)	
	# Delete all of the old generations information.
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	hindex = np.delete(np.asarray(hindex,dtype=float),deleteallindex)
	# Keep genes around if within burnin faze
	if gen >= geneswap:
		genes = np.delete(genes,deleteallindex,axis=0)
		genes = np.asarray(genes,dtype=int)
		genes = genes.tolist()
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)
	
	# Just one more shuffle to mix up empty versus killed off grids
	shuffle(freegrid)

	# Store total number of Deaths information for output 
	# Check that agedeaths equals age distribution file
	Track_MDeaths.append((np.asarray(agedeaths)/2.).tolist())
	Track_FDeaths.append((np.asarray(agedeaths)/2.).tolist())
	if len(extra_agedeaths) != 0:
		Track_MDeaths[gen][-1] = Track_MDeaths[gen][-1] + sum(extra_agedeaths)/2.
		Track_FDeaths[gen][-1] = Track_FDeaths[gen][-1] + sum(extra_agedeaths)/2.
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,mature,hindex
	return tupMort
	# End::DDMortality
	
# ---------------------------------------------------------------------------------------------------	 
def AdultSelection(tupMort,fitvals,mature,SelectionDeaths):
	'''
	AdultSelection()
	Mortality of old generation
	using selection values. Only for mature individuals
	'''
	# Unpack tuple
	freegrid = tupMort[0]
	id = tupMort[1]
	sex = tupMort[2]
	age = tupMort[3]
	xgrid = tupMort[4]
	ygrid = tupMort[5]
	genes = tupMort[6]	
	FID = tupMort[7]
	infection = tupMort[8]
	mature = tupMort[9]
	hindex = tupMort[10]
	
	deleteallindex = []
	# Loop through mature individuals
	for i in range(len(mature)):
		if mature[i] == 1:
			# Find it's location
			usefitvals = fitvals[FID[i]]
			# Check genotype and match fitvals
			if int(genes[i][0][0]) == 2: # AA
				diffmort = float(usefitvals[0])/100.
			elif int(genes[i][0][0]) == 1 and int(genes[i][0][1]) == 1: # Aa
				diffmort = float(usefitvals[1])/100.
			elif int(genes[i][0][1]) == 2: #aa
				diffmort = float(usefitvals[2])/100.
			else: # Another genotype
				diffmort = 0.0
			
			# Then flip the coin to see if  individual survives its location
			randcheck = np.random.uniform()			
			# If individual did not survive: keep track of delete ones
			if randcheck < diffmort:
				SelectionDeaths.append(1) # Record
				# Keep storage spots to delete
				deleteallindex.append(i)				
				# Add to freegrid
				freegrid.append(FID[i])
	
	# Delete spots
	deleteallindex = np.asarray(deleteallindex,dtype='int')			
	sex = np.delete(sex,deleteallindex)
	id = np.delete(id,deleteallindex)
	age = np.delete(age,deleteallindex)
	xgrid = np.delete(xgrid,deleteallindex)
	ygrid = np.delete(ygrid,deleteallindex)
	mature = np.delete(mature,deleteallindex)
	genes = np.delete(genes,deleteallindex,axis=0)
	FID = np.delete(FID,deleteallindex)
	infection = np.delete(np.asarray(infection),deleteallindex)	
	hindex = np.delete(np.asarray(hindex,dtype=float),deleteallindex)		
	
	# Return variables from this argument
	tupMort = freegrid,id,sex,age,xgrid,ygrid,genes,FID,infection,hindex
	return tupMort
	# End::AdultSelection
	
# ---------------------------------------------------------------------------------------------------	 
def DoMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,FID,Magemort,Fagemort,infection,geneswap,popmodel,K_env,fitvals,mature,cdevolveans,Opt3SelectionDeaths,startSelection,subpop,hindex):
	'''
	DoMortality()
	Mortality of old generation
	Input: Adult mortality% 
	Output: Old files minus the killed off individuals:
	freegrid = [xycdmatid location of free grid spot 
	in random order]
	Constant or Density dependent functions here.
	'''
	
	Opt3SelectionDeaths.append([]) # Spot added for generation
	# Switch for model choice
	if popmodel == 'exp':
	
		tupMort = ConstantMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,FID,Magemort,Fagemort,infection,geneswap,mature,subpop,hindex)
		
	elif popmodel == 'logistic' or popmodel == 'rickers' or popmodel == 'richards':
		
		tupMort = DDMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,FID,infection,geneswap,K_env,popmodel,Magemort,Fagemort,mature,subpop,hindex)	
		
	else:
		print('This population model for survival does not exist')
		sys.exit(-1)
	
	# Then apply fitness values to adults if there and mature
	if cdevolveans == '3' and gen >= startSelection:
		tupMort = AdultSelection(tupMort,fitvals,mature,Opt3SelectionDeaths[gen])
	Opt3SelectionDeaths[gen] = sum(Opt3SelectionDeaths[gen]) 
		
	return tupMort
	
	# End::DoAdultMortality()
	
# ---------------------------------------------------------------------------------------------------	 
def AddIndividuals(cdclimgentime,gen,id_here,age_here,genes_here,sex_here,FID_here,infection_here,allelst,xyfilename,datadir,alleles,hindex_here,subpop,freegrid):
	'''
	AddIndividuals()
	This function reads in the multiple xy files given at specified cdclimate year and adds individuals to each subpopulation. Checks for more individuals than K are given. 
	'''
	
	# This is the index into the xyfilename and allefilename 
	ThisIndex = np.where(np.asarray(cdclimgentime) == str(gen))[0][0]
	
	# Read in XY File
	xy = ReadXY(datadir+xyfilename[ThisIndex])
	
	# Just a quick check, should be 19 long
	if len(xy[0]) != 20:
		print('XY files should be 20 long.')
		sys.exit(-1)
	
	# Store information (some not needed, e.g., xgrid, ygrid). 
	FID_add = [] # store index of each added individual
	subpop_add = [] # store subpop id of each added individual
	id_add = []
	sex_add = []	
	age_add = []	
	genes_add = [] # For known genes, need to update below
	infection_add = []
	# gridmort is read in through cdclimate in case changing there; first file will then give values through time; update this?
	
	# Read in new added individuals information here
	for i in range(len(xy)-1):
		FID_add.append(str(i))
		subpop_add.append(xy[i+1][0])
		id_add.append(xy[i+1][4])
		sex_add.append(xy[i+1][5])
		age_add.append(xy[i+1][6])
		infection_add.append(xy[i+1][7])		
		# Error check for age in known file must be 1 or greater
		if age_add[i] != 'NA':
			if int(age_add[i]) < 1:
				print('Known file must initize with age 1+.')
				sys.exit(-1)
		# Change id, subpop to 'NA', FID to 'NA'
		if sex_add[i] == 'NA':
			subpop_add[i] = 'NA'
			FID_add[i] = 'NA'
		
	# turn new and add into arrays
	subpop_here = np.asarray(subpop)[FID_here] # Since subpop doesn't change, then extract IDs for each location here
	freegrid = np.asarray(freegrid)
	
	# Get K for each pop or location (1)
	subpopK = count_unique(np.asarray(subpop))	
	
	# Add the new individuals to each pop (or location/grid spot), check if space
	for isub in range(len(subpopK[0])):
		
		thisPop = subpopK[0][isub] # make sure on right pop given strings
		thisK = subpopK[1][isub] # grab K for this pop / or location/grid
		
		# How many individuals are currently in this subpop or lcoation/grid
		#count_currentN = len(np.where(sex_here[np.where(subpop_here == thisPop)[0]] != 'NA')[0])
		#count_currentN = len(np.where(sex_here[np.where(FID_here == int(thisPop)-1)[0]] != 'NA')[0]) # Not sure works with pops now
		
		# Count how many individuals in this spot
		count_currentN = len(np.where(subpop_here == thisPop)[0])
		
		# How many individuals are going to be added to this pop or spot
		addN = np.where(np.asarray(subpop_add) == thisPop)[0] # could be > 1, index / FID location of added individuals
		count_addN = len(addN)
		
		# If there are no individuals to add, then continue
		if count_addN == 0:
			continue
		else: # Else keep going and adding individuals		
			# Check to see if over K?
			if thisK < count_currentN + count_addN:
				print(('Warning: Exceeded carrying capacity when adding individuals to this subpopulation '+thisPop+ ' and no individual added here for this gen '+str(gen))) 
			else:
				# Loop through each ind to add; could be > 1
				for iadd in range(count_addN):
					thisAddInd = addN[iadd] # index/FID of added individual
					# Make sure there is not an individual already here
					checkAddInd = np.where(thisAddInd == np.asarray(FID_here))[0]
					if len(checkAddInd) > 0: # There is ind here, move onto next add individual in this pop, else move to next subpop
						#pdb.set_trace()
						continue
					else: # No one here, continue to adding andupdate the *_here variables
						id_here.append(id_add[thisAddInd])
						sex_here.append(sex_add[thisAddInd])
						age_here.append(age_add[thisAddInd])
						infection_here.append(infection_add[thisAddInd])
						FID_here.append(int(FID_add[thisAddInd]))
						freegrid = np.delete(freegrid,np.where(thisAddInd==freegrid)[0])	# Delete from freegrid					
													
						# Then update the genes_here spots with allele freq file							
						genes_here.append([])	
						# For each loci:
						for j in range(len(allelst[ThisIndex])):							
							# Take a random draw from the w_choice function at jth locus
							# Using the first allele file in list
							rand1 = w_choice_general(allelst[ThisIndex][j])[0]
							rand2 = w_choice_general(allelst[ThisIndex][j])[0]
							
							# Append assinment onto indall array - run through each condition for assignment of 1s or 2s or 0s
							# 	1s = heterozygous at that locus
							#	2s = homozygous at that locus
							#	0s = absence of allele
							for k in range(len(allelst[ThisIndex][j])):									
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
								genes_here[-1].append(str(tempindall))	
								
						# ---------------------------------------------
						# Update Hindex
						# ---------------------------------------------
						if genes_here[-1][0] == str(2):
							hindex_here.append(str(1.0))
						elif genes_here[-1][1] == str(2):
							hindex_here.append(str(0.0))
						elif genes_here[-1][0] == 1 and genes_here[-1][1] == 1:
							hindex_here.append(str(0.5))
						else:
							hindex_here.append(str(-9999))	
				
	# Turn back into lists?	
	return id_here,sex_here,age_here,genes_here,infection_here,hindex_here,FID_here,freegrid.tolist()
	
	# End::AddIndividuals()

# ---------------------------------------------------------------------------------------------------	 
def DoEpigenetics(epimod,betas,sex,id,age,genes,infection,Track_EpigeneMod1,Track_EpigeneMod2,Track_EpigeneDeaths,gen,cdevolveans,epigeneans,startEpigene,geneswap,alleles,loci):	
	'''
	The function modifies the genotype region specified with given probability at location individual
	is at. First loci are fixed DNA changes with Selection module, then next loci are epigenetic regions, the rest are neutral.
	'''
	
	# Tracking
	Track_EpigeneMod1.append([])
	Track_EpigeneMod2.append([])
	Track_EpigeneDeaths.append([])
	
	# Get location in genes array for epigene region
	# ----------------------------------------------
	genes = np.asarray(genes)
	if cdevolveans.split('_')[0] == 'M':
		# Then the first l loci are for selection, and next epigenetics
		selloci = int(cdevolveans.split('_')[2].split('L')[1])
		epiloci = int(epigeneans.split('_')[1].split('L')[1])
	else:
		selloci = 0
		epiloci = int(epigeneans.split('_')[1].split('L')[1])
	epiloci_index = list(range(selloci*int(cdevolveans.split('_')[3].split('A')[1]),selloci*int(cdevolveans.split('_')[3].split('A')[1])+epiloci*int(epigeneans.split('_')[2].split('A')[1]))) # Index location for epi region, epi region only considers 2 alleles but generic form to index to here.	
	 
	# If first generation, set genes epigenetic region to 0 or turned off, unless geneswap not started
	if gen >= geneswap:
		genes[:,epiloci_index] = 0
	
	# Skip if delayed start time
	if gen >= startEpigene:
		deleteallindex = []
		
		# Loop through individuals - get epimutation
		# ------------------------------------------
		for iind in range(len(epimod)):	
			# Skip if no one at this spot
			if sex[iind] != 'NA':
				linmodel = [] # [loci][alleles]	for linear model
				# Loop through each locus position
				#for ilocus in xrange(len(epiloci_index)):
				for ilocus in range(epiloci):
					
					#thisregion = genes[iind][epiloci_index][(ilocus*int(epigeneans.split('_')[2].split('A')[1])):(ilocus*int(epigeneans.split('_')[2].split('A')[1])+2)]
					
					
					# If alleles are [0,2] or [2,0] then skip this step
					if int(genes[iind][epiloci_index][ilocus*2]) == 2 or int(genes[iind][epiloci_index][ilocus*2+1]) == 2:
						continue
					
					else:
						# Check first allele
						# -----------------------------------
						if int(genes[iind][epiloci_index][ilocus*2]) == 0: # Only if first allele off
							# Check for modification here
							epimutateprob = float(epimod[iind][ilocus].split(';')[0])
							randno = np.random.uniform()
							if randno < epimutateprob: # mutation occurs, turn on
								get_epialleleindex = epiloci_index[ilocus*2]
								genes[iind][get_epialleleindex] = 1
								Track_EpigeneMod1[gen].append(1)
							else: # No mutation
								Track_EpigeneMod1[gen].append(0)
						
						# Then check for dependence (second allele)
						# ----------------------------------
						if epigeneans.split('_')[4] == 'Dep':
							# Then if first allele is on or just got turned on
							if int(genes[iind][epiloci_index][ilocus*2]) != 0:
								# And if the second allele is not on
								if int(genes[iind][epiloci_index][ilcous*2+1]) == 0:	
									# Check for methylation modifciation here
									epimethylmutateprob = float(epimod[iind][ilocus].split(';')[1])
									randno = np.random.uniform()
									if randno < epimethylmutateprob: # mutation occurs
										get_epialleleindex = epiloci_index[ilocus*2+1]
										genes[iind][get_epialleleindex] = 1
										Track_EpigeneMod2[gen].append(1)
									else: # No mutation
										Track_EpigeneMod2[gen].append(0)
								else: 
									Track_EpigeneMod2[gen].append(0)
										
						# Or Independence (second allele)
						# -------------------------------
						elif epigeneans.split('_')[4] == 'Ind':
							# For the second allele, if it is off
							if int(genes[iind][epiloci_index][ilocus*2+1]) == 0:
								# Check for modification here
								epimutateprob = float(epimod[iind][ilocus].split(';')[1])
								randno = np.random.uniform()
								if randno < epimutateprob: # mutation occurs, turn on
									get_epialleleindex = epiloci_index[ilocus*2+1]
									genes[iind][get_epialleleindex] = 1
									Track_EpigeneMod2[gen].append(1)
								else: # No mutation
									Track_EpigeneMod2[gen].append(0)
													
						# OR something entered wrong, check
						else:
							print('Epigenetic answer needs to specify Ind or Dep, see usermanual examples.')
							sys.exit(-1)
					'''
					# Next Apply fitness consequence for this locus
					# ---------------------------------------------
					if epigeneans.split('_')[3] == 'ModelY': # Code 1,0
						if int(genes[iind][epiloci_index[ilocus]][0]) == 2:
							allele1 = 1
						else:
							allele1 = int(genes[iind][epiloci_index[ilocus]][0])
						if int(genes[iind][epiloci_index[ilocus]][1]) == 2:
							allele2 = 1
						else:
							allele2 = int(genes[iind][epiloci_index[ilocus]][1])
					elif epigeneans.split('_')[3] == 'ModelX': # Code 2,1,0
						allele1 = int(genes[iind][epiloci_index[ilocus]][0])
						allele2 = int(genes[iind][epiloci_index[ilocus]][1])
					
					linmodel.append(betas[ilocus][0]*allele1 + betas[ilocus][1]*allele2)	
					'''		
				'''
				# Add the beta not to linear model after loci loop
				linmodel.append(betas[-1])
		
				# Add the linear model together and logit
				#Fitness = np.exp(sum(linmodel)) / (1. + np.exp(sum(linmodel)))

				# Convert fitness to differential mortality - 1 - Fitness
				diffmort = 1. - Fitness
				
				# Check if this individual survives
				randno = np.random.uniform()
				if randno < diffmort: # Did not survive
					Track_EpigeneDeaths[gen].append(1)
					# Keep sotrage spots to delete
					deleteallindex.append(iind)
				'''
		
		# Add NA to spots (sex,age,genes,infection), id gets 'OPEN'
		#deleteallindex = np.asarray(deleteallindex,dtype='int')
		id = np.asarray(id)
		#id[deleteallindex] = 'OPEN'
		sex = np.asarray(sex,dtype='|U2')
		#sex[deleteallindex] = 'NA'
		age = np.asarray(age,dtype='|U2')
		#age[deleteallindex] = 'NA'
		infection = np.asarray(infection,dtype='|U2')
		#infection[deleteallindex] = 'NA'
		genes = np.asarray(genes,dtype='|U2')
		#genes[deleteallindex] = 'NA'
			
		# Summary Tracking numbers, cleanup
		Track_EpigeneMod1[gen] = sum(Track_EpigeneMod1[gen])
		Track_EpigeneMod2[gen] = sum(Track_EpigeneMod2[gen])
		Track_EpigeneDeaths[gen] = sum(Track_EpigeneDeaths[gen])
		return id.tolist(),sex.tolist(),age.tolist(),genes.tolist(),infection.tolist()
	else: # return original lists
		Track_EpigeneMod1[gen] = sum(Track_EpigeneMod1[gen])
		Track_EpigeneMod2[gen] = sum(Track_EpigeneMod2[gen])
		Track_EpigeneDeaths[gen] = sum(Track_EpigeneDeaths[gen])
		return id,sex,age,genes.tolist(),infection
	 
	# End::DoEpigenetics()