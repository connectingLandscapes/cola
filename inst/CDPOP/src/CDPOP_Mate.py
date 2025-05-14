# -------------------------------------------------------------------------------------------------
# CDPOP_Mate.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for mate processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError("Numpy required.")
	
# Python specific functions
import pdb, random, os, sys, copy
#from sets import Set

# --------------------------------------------------------------------------
def countDuplicatesInList(dupedList):
	'''
	countDuplicatesInList() - Counts dupicates in lists
	'''
	uniqueSet = set(item for item in dupedList)
	return [dupedList.count(item) for item in uniqueSet]
	
	# End::countDuplicatesInList()

# ---------------------------------------------------------------------------------------------------
def count_unique(keys):
    uniq_keys = np.unique(keys)
    bins = uniq_keys.searchsorted(keys)
    return uniq_keys, np.bincount(bins)
	
	#End::count_unique()
	
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
def DoSexualNY(nomales,xycdmatrix,females,count,\
males,matemovethresh,Bearpairs,subpop,gen,selfans,matefreq):
	'''
	DoSexualNY()
	This function is the mating function for:
	sexual reproduction
	females	without replacement
	males with replacement
	'''	
	
	# Create an empty probability array to be appended to
	probarray = []
	
	# Check if this female will mate
	randmate = np.random.uniform()
	
	#Mating occurs or potential to occur
	if randmate < matefreq:
						
		# Make array of individuals, removing itself unless selfing on
		if selfans == 'Y':
			indspots = np.asarray(copy.deepcopy(males))
		else:
			indspots = np.asarray(copy.deepcopy(males))
			delspot = np.where(indspots==females[count])[0]
			indspots = np.delete(indspots,delspot)
		shuffle(indspots)
		
		# Extract each male grid probability
		probarray.append(xycdmatrix[females[count]][indspots])
				
		# If statement to check if there were individuals in probarray:
		if sum(probarray[0]) != 0.0:
		
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
		
			# And store the mated pair information.						
			Bearpairs.append([females[count],indspots[itemselect]])
								
		# If statement to check if there were not males < matemovethresh:
		elif sum(probarray[0]) == 0.0:
			
			# Store the female that did not mate with -9999 designation
			Bearpairs.append([females[count],-9999])

	# Mating does not happen this year
	else:		
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[count],-9999])
	
	return Bearpairs,males
	
	# End::DoSexualNY()
	
# ---------------------------------------------------------------------------------------------------	
def DoSexualYY(nomales,xycdmatrix,females,males,\
matemovethresh,Bearpairs,nofemales,subpop,gen,selfans,matefreq):
	'''
	DoSexualYY()
	This function is the mating function for: 
	sexual and asexual reproduction
	females	with replacement
	males with replacement.
	'''	
	# Create an empty probability array to be appended to
	probarray = []					

	# Randomly grab a female
	intfemale = int((nofemales)*np.random.uniform())
	
	# Check if this female will mate
	randmate = np.random.uniform()
	
	#Mating occurs or potential to occur
	if randmate < matefreq:
	
		# Make array of individuals, removing itself unless selfing on
		if selfans == 'Y':
			indspots = np.asarray(copy.deepcopy(males))
		else:
			indspots = np.asarray(copy.deepcopy(males))
			delspot = np.where(indspots==females[intfemale])[0]
			indspots = np.delete(indspots,delspot)
		shuffle(indspots)
		
		# Extract each male grid probability
		probarray.append(xycdmatrix[females[intfemale]][indspots])
		
		# If statement to check if there were individuals in probarray:
		if sum(probarray[0]) != 0.0:
			
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
				
			# And store the mated pair information.						
			Bearpairs.append([females[intfemale],indspots[itemselect]])
							
		# If statement to check if there were not males < matemovethresh:
		elif sum(probarray[0]) == 0.0:
		
			# Store the female that did not mate with -9999 designation
			Bearpairs.append([females[intfemale],-9999])				
	
	# Mating does not happen this year
	else:		
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[intfemale],-9999])
	
	# Return Variables from this function
	return Bearpairs,males
	
	# End::DoSexualYY()

# ---------------------------------------------------------------------------------------------------	
def DoSexualNN(nomales,xycdmatrix,females,count,\
males,matemovethresh,Bearpairs,subpop,gen,matefreq):
	'''
	DoSexualNN()
	This function is the mating function for
	sexual reproduction
	females	with replacement
	males with replacement
	'''	
	# Create an empty probability array to be appended to
	probarray = []
					
	# Check if this female will mate
	randmate = np.random.uniform()
	
	#Mating occurs or potential to occur
	if randmate < matefreq:
	
		# Extract each male grid probability
		probarray.append(xycdmatrix[females[count]][males])
					
		# If statement to check if there were individuals in probarray:
		if sum(probarray[0]) != 0.0:
			
			# Select the w_choice item
			itemselect = w_choice_item(probarray[0])
		
			# And store the mated pair information.						
			Bearpairs.append([females[count],males[itemselect]])
			
			# Then delete that male from the male list
			males = np.delete(males,itemselect)
								
		# If statement to check if there were not males < matemovethresh:
		elif sum(probarray[0]) == 0.0:
		
			# Store the female that did not mate with -9999 designation
			Bearpairs.append([females[count],-9999])	
			
	# Mating does not happen this year
	else:		
		# Store the female that did not mate with -9999 designation
		Bearpairs.append([females[count],-9999])
	
	return Bearpairs,males
	
	# End::DoSexualNN()	

# ---------------------------------------------------------------------------------------------------	 
def DoMate(nogrids,sex,age,\
freplace,mreplace,matemoveno,matemovethresh,\
xycdmatrix,MateDistED,\
MateDistCD,xgridcopy,ygridcopy,ToTMales,\
ToTFemales,BreedMales,BreedFemales,\
sexans,selfans,\
MateDistEDstd,MateDistCDstd,FAvgMate,MAvgMate,\
FSDMate,MSDMate,filledgrids,Female_BreedEvents,\
gen,subpop,BreedFemales_age,agemort,Mmature,Fmature,ScaleMax,ScaleMin,A,B,C,MateDistances,matefreq):

	'''
	DoMate()
	This is the mating function for choosing
	individual mate pairs. Two calls here:
	sexual and asexual mating.	
	'''
	
	# --------------------------------------------------------
	# Preliminary: Needed for both sexual and asexual routines	
	# --------------------------------------------------------
	
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	unique_subpops = np.unique(subpop)
	mature = []
	
	# ---------------------------------------
	# Step10a: Call DoMateSexual()
	# ---------------------------------------
	if (sexans == 'Y'):
		'''
		DoMateSexual()
		This function is the mating function for 
		sexual reproduction.
		'''
		# ---------------------------------------------------
		# Select males and females for mating
		# ---------------------------------------------------
		
		# Check for mature individuals
		# ----------------------------
		allfemales = []
		allmales = []
		females = [] # Breeding age females
		males = [] # Breeding age males
		for i in range(nogrids):
			Indsex = sex[i]
			Indage = age[i]
			# If there are multiple subpop maturation rates given
			if len(Fmature) > 1:
				Indpop = int(subpop[i])-1 # Subtract 1 for count 0
			else:
				Indpop = 0
			if Indsex != 'NA':
				# Last age check
				if int(Indage) >= len(Fmature[Indpop]):
					Indage = str(len(Fmature[Indpop]) - 1)
				# If Female
				if str(Indsex) == '0': 
					allfemales.append(i)
					# Carful of indexing Fmature 1 +
					matval = Fmature[Indpop][int(Indage)-1]
					randmat = np.random.uniform()
					# Becomes mature
					if randmat < matval: 
						females.append(i)
						mature.append(1)
					else:
						mature.append(0)
				# If male
				elif str(Indsex) == '1':
					allmales.append(i)
					matval = Mmature[Indpop][int(Indage)-1]
					randmat = np.random.uniform()
					# Becomes mature
					if randmat < matval: 
						males.append(i)
						mature.append(1)
					else:
						mature.append(0)
			else:
				mature.append('NA')
		# Then get the length of each sex that are reproducing
		nomales = len(males)
		nofemales = len(females)
		
		# And then sum them up - and store numbers
		ToTMales.append([])
		ToTFemales.append([])
		BreedMales.append([])
		BreedFemales.append([])
		BreedFemales_age.append([])
		ToTMales[gen].append(len(allmales))
		ToTFemales[gen].append(len(allfemales))
		BreedMales[gen].append(nomales)
		BreedFemales[gen].append(nofemales)
		for i in range(nosubpops):
			# Get subpop location 
			popall = np.where(np.asarray(subpop) == str(i+1))[0]
			# Males in pop
			popmales = np.in1d(popall,allmales)
			ToTMales[gen].append(sum(popmales))
			# Females in pop
			popfemales = np.in1d(popall,allfemales)
			ToTFemales[gen].append(sum(popfemales))
			# Breed males in pop
			popmales = np.in1d(popall,males)
			BreedMales[gen].append(sum(popmales))
			# Breed females in pop
			popfemales = np.in1d(popall,females)
			BreedFemales[gen].append(sum(popfemales))
			
		# Get count of female age breeders
		for i in range(len(agemort[0])):
			BreedFemales_age[gen].append([])
				
		# Check if females left
		if len(females) != 0:
			Fsexreproage = np.asarray(age)[np.asarray(females)]
			countFage = count_unique(Fsexreproage)		
			count = 0
			for item in countFage[0]:
				if item != 'NA':
					if int(item) > len(agemort[0]):
						item = len(agemort[0])-1
						BreedFemales_age[gen][int(item)-1].append(countFage[1][count])
					else:
						BreedFemales_age[gen][int(item)-1].append(countFage[1][count])
				count = count+1
		# Sum up population age tracker
		for i in range(len(agemort[0])):
			BreedFemales_age[gen][i] = sum(BreedFemales_age[gen][i])
		
		# Choose mate for each female or individual
		Bearpairs = []	# Empty matrix: xy indexes
		
		# Shuffle the females
		shuffle(females)
		
		# If there were no reproducing males or females
		if nomales == 0 or nofemales == 0:
			Bearpairs.append([-9999,-9999])
		
		# If there were reproducing males and females
		if nomales != 0 and nofemales != 0:
		
			# For the case of a Female without replacement and a male with replacement
			if freplace == 'N' and mreplace == 'Y':
							
				# Loop through while loop until all females paired up.		
				count = 0		# Initialize the while loop
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)
				while count < nofemales:
								
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNY(nomales,xycdmatrix,females,count,\
					tempmales,matemovethresh,Bearpairs,subpop,gen,selfans,matefreq)
						
					# Update count
					count = count + 1
					
			# For the case of a Female with replacement and a male with replacement
			elif freplace == 'Y' and mreplace == 'Y':
					
				# Loop through while loop until all females paired up, but do this nogrid times.		
				count = 0		# Initialize the while loop
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)
				while count < filledgrids:
						
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualYY(nomales,xycdmatrix,females,\
					tempmales,matemovethresh,Bearpairs,nofemales,subpop,gen,selfans,matefreq)
														
					# Update count
					count = count + 1
									
			# For the case of Female with replacement and male without replacement
			elif freplace == 'Y' and mreplace == 'N':
			
				print('Female with replacement and Male without replacement not coded yet.')
				sys.exit(-1)
				
			# For the case of Female without replacement and male without replacement
			elif freplace == 'N' and mreplace == 'N':
				
				# Loop through while loop until all male female pairs occur		
				count = 0		# Initialize the while loop for females
				# Create a temp male to delete from
				tempmales = copy.deepcopy(males)				
				while count < nofemales:
								
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNN(nomales,xycdmatrix,females,count,tempmales,matemovethresh,Bearpairs,subpop,gen,matefreq)
											
					# Update count
					count = count + 1
			
			# Error check
			else:
				print('This Female/Male mating structure does not exist. Must be Y/N combinations.')
				sys.exit(-1)		
		
		# Sort the females for later functions
		females = list(np.sort(females))
		
	# ---------------------------------------
	# Step10b: Call DoMateAsexual()
	# ---------------------------------------
	if (sexans=='N'):	
		'''
		DoMateAsexual()
		This function is the mating function for 
		asexual reproduction.
		'''
			
		# ---------------------------------------------------
		# Select individuals for mating
		# ---------------------------------------------------
		
		# Check for mature individuals
		# ----------------------------
		allgrids = []
		breedgrids = []	
		mature = []
		# Loop through and grab each index for each catagory for reproage
		for i in range(nogrids):
			Indsex = sex[i]
			Indage = age[i]
			# If there are multiple subpop maturation rates given
			if len(Fmature) > 1:
				Indpop = int(subpop[i])-1 # Subtract 1 for count 0
			else:
				Indpop = 0
			if Indsex != 'NA':
				# Last age check
				if int(Indage) >= len(Fmature[Indpop]):
					Indage = str(len(Fmature[Indpop]) - 1)
				# If Female (although doesn't matter, but looks up Fmature)
				if str(Indsex) == '0': 
					allgrids.append(i)
					# Carful of indexing Fmature 1 +
					matval = Fmature[Indpop][int(Indage)-1]
					randmat = np.random.uniform()
					# Becomes mature
					if randmat < matval: 
						breedgrids.append(i)
						mature.append(1)
					else:
						mature.append(0)
				# If male (again doesn't matter, but looks up Mmature)
				elif str(Indsex) == '1':
					allgrids.append(i)
					matval = Mmature[Indpop][int(Indage)-1]
					randmat = np.random.uniform()
					# Becomes mature
					if randmat < matval: 
						breedgrids.append(i)
						mature.append(1)
					else:
						mature.append(0)
			else:
				mature.append('NA')
							
		# Then get the length of each sex that are reproducing
		nobreedgrids = len(breedgrids)
		
		# And then sum them up - and store numbers
		ToTMales.append([])
		ToTFemales.append([])
		BreedMales.append([])
		BreedFemales.append([])
		BreedFemales_age.append([])
		ToTMales[gen].append(len(allgrids))
		ToTFemales[gen].append(len(allgrids))
		BreedMales[gen].append(nobreedgrids)
		BreedFemales[gen].append(nobreedgrids)
		for i in range(nosubpops):
			# Get subpop location 
			popall = np.where(np.asarray(subpop) == str(i+1))[0]
			# Males in pop
			popmales = np.in1d(popall,allgrids)
			ToTMales[gen].append(sum(popmales))
			# Females in pop
			popfemales = np.in1d(popall,allgrids)
			ToTFemales[gen].append(sum(popfemales))
			# Breed males in pop
			popmales = np.in1d(popall,breedgrids)
			BreedMales[gen].append(sum(popmales))
			# Breed females in pop
			popfemales = np.in1d(popall,breedgrids)
			BreedFemales[gen].append(sum(popfemales))
		
		# Get count of breeders
		for i in range(len(agemort[0])):
			BreedFemales_age[gen].append([])
		
		# Check if breeders left
		if len(breedgrids) != 0:
			Fsexreproage = np.asarray(age)[np.asarray(breedgrids)]
			countFage = count_unique(Fsexreproage)		
			count = 0
			for item in countFage[0]:
				if item != 'NA':
					if int(item) > len(agemort[0]):
						item = len(agemort[0])-1
						BreedFemales_age[gen][int(item)-1].append(countFage[1][count])
					else:
						BreedFemales_age[gen][int(item)-1].append(countFage[1][count])
				count = count+1
		# Sum up population age tracker
		for i in range(len(agemort[0])):
			BreedFemales_age[gen][i] = sum(BreedFemales_age[gen][i])
				
		# Choose mate for each female or individual
		Bearpairs = []	# Empty matrix: xy indexes
		
		# Shuffle the grids
		shuffle(breedgrids)
		
		# If there were no reproducing grids
		if nobreedgrids == 0:
			Bearpairs.append([-9999,-9999])	

		# If there were reproducing grids
		else:
					
			# For the case of a first mate without replacement and second mate with replacement
			if freplace == 'N' and mreplace == 'Y':
			
				# Loop through while loop until all grids are paired up according to function		
				count = 0		# Initialize the while loop
				# Create a temp grid to delete from
				tempmales = copy.deepcopy(breedgrids)
				while count < nobreedgrids:
													
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualNY(nobreedgrids,xycdmatrix,breedgrids,count,\
					tempmales,matemovethresh,Bearpairs,subpop,gen,selfans,matefreq)
																			
					# Update count
					count = count + 1
					
			# For the case of a first mate with replacement and second mate with replacement
			elif freplace == 'Y' and mreplace == 'Y':
			
				# Loop through filledgrids times.		
				count = 0		# Initialize the while loop
				# Create a temp male to delete from
				tempmales = copy.deepcopy(breedgrids)
				while count < filledgrids:
				
					# Get probability function of user defined input number
					Bearpairs,tempmales = DoSexualYY(nobreedgrids,xycdmatrix,breedgrids,\
					tempmales,matemovethresh,Bearpairs,nobreedgrids,subpop,gen,selfans,matefreq)
										
					# Update count
					count = count + 1					
				
			# For the case of first mate without replacement and second mate without replacement
			elif freplace == 'N' and mreplace == 'N':
			
				print('This mate function is not coded up yet, email Erin.')
				sys.exit(-1)
				
			# For the case of the first mate with replacement and second mate without
			elif freplace == 'Y' and mreplace == 'N':				
				print('This is not coded up yet, but is the same case as NY replacement.')
				sys.exit(-1)
				
			# Error check
			else:
				print('This asexual mating structure does not exist. Must be Y/N combinations.')
				sys.exit(-1)
		
		# Sort the females for later functions
		females = list(np.sort(breedgrids))
		males = list(np.sort(breedgrids))
		
	# ----------------------------------------
	# Summary Stats on Mate functions
	# ----------------------------------------
	# Store the average distance mates were choosen 
	# temp variable to store the number of -9999 numbers
	tempImmStorage = []
	# temp variable to store mate distance
	tempAvgMateED = []
	tempAvgMateCD = []
		
	# Loop through each CDpair
	for ipair in range(len(Bearpairs)):
		
		# if -9999 store
		if Bearpairs[ipair][1]==-9999:
			tempImmStorage.append(1)
		
		# else calculate average mate distance
		else:
			tempAvgMateED.append(np.sqrt((xgridcopy[Bearpairs[ipair][0]]-xgridcopy[Bearpairs[ipair][1]])**2+(ygridcopy[Bearpairs[ipair][0]]-ygridcopy[Bearpairs[ipair][1]])**2))
			probval = xycdmatrix[Bearpairs[ipair][0]][Bearpairs[ipair][1]]
			
			if matemoveno == '4' or matemoveno == '6':
				cdval = 0.0
				
			# If linear
			elif matemoveno == '1':
				cdval = (probval - 1.) * (-matemovethresh) 
				
			# if inverse square
			elif matemoveno == '2':
				probval = xycdmatrix[Bearpairs[ipair][0]][Bearpairs[ipair][1]]
				if probval == 1.0:
					cdval = 0.0
				else:
					cdval = np.sqrt(1. / (probval * (ScaleMax - ScaleMin) + ScaleMin))
				
			# If neg exp
			elif matemoveno == '5':
				cdval = np.log((probval * (ScaleMax-ScaleMin) + ScaleMin)/float(A)) / (-float(B) * np.log(10))
			
			# If gaussian
			elif matemoveno == '7':				
				cdval = float(B) + np.sqrt(-2*float(C)**2 * np.log((probval*(ScaleMax-ScaleMin)+ScaleMin)/float(A)))
			
			# If just rescaled
			elif matemoveno == '8':
				cdval = probval*(ScaleMax-ScaleMin)+ScaleMin

			# If straight prob matrix
			elif matemoveno == '9':
				cdval = probval
				
			else:
				print('Mate move function does not exist.')
				sys.exit(-1)
			tempAvgMateCD.append(cdval)		
	
	# If at least some individuals mated
	if len(Bearpairs) > sum(tempImmStorage):
		
		# And append to MateDistED
		MateDistED.append(sum(tempAvgMateED)/(len(Bearpairs)-sum(tempImmStorage)))
		MateDistEDstd.append(np.std(tempAvgMateED))

		# And append to MateDistCD
		MateDistCD.append(sum(tempAvgMateCD)/(len(Bearpairs)-sum(tempImmStorage)))
		MateDistCDstd.append(np.std(tempAvgMateCD))
		
	# If all immigrants
	else:
	
		# And append to MateDistED
		MateDistED.append(0)
		MateDistEDstd.append(0)

		# And append to MateDistCD
		MateDistCD.append(0)
		MateDistCDstd.append(0) 
	
	# Track actual number of breeding events of females.
	Female_BreedEvents.append(len(Bearpairs)-sum(tempImmStorage))
	MateDistances.append(tempAvgMateCD)	
	
	# Delete temp storage
	del(tempImmStorage)
	del(tempAvgMateED)
	del(tempAvgMateCD)
		
	# Get the avgerage and sdev for female and male mating times
	#	For variance in reproductive succes.
	Bearpairs_array = np.asarray(Bearpairs)
	femalesmated = Bearpairs_array[:,0][np.where(Bearpairs_array[:,1]>=0)]
	femalesmated = list(femalesmated)
	malesmated = Bearpairs_array[:,1][np.where(Bearpairs_array[:,1]>=0)]
	malesmated = list(malesmated)
	females_nomate = list(Bearpairs_array[:,0][np.where(Bearpairs_array[:,1]==-9999)])
	males_nomate = list(Bearpairs_array[:,1][np.where(Bearpairs_array[:,0]==-9999)])
	# Append to tracker
	FAvgMate.append(np.mean(countDuplicatesInList(femalesmated)))
	MAvgMate.append(np.mean(countDuplicatesInList(malesmated)))
	FSDMate.append(np.std(countDuplicatesInList(femalesmated)))
	MSDMate.append(np.std(countDuplicatesInList(malesmated)))
	
	# Return variables from this function
	tupMate = Bearpairs, females, females_nomate, males, males_nomate,mature
	return tupMate
	
	#End::DoMate()