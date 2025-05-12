# -------------------------------------------------------------------------------------------------
# CDPOP_Offspring.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file for offspring processes.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	from numpy.random import *
	import numpy as np
except ImportError:
	raise ImportError("Numpy required.")
import pdb,sys
#from sets import Set
import random
from collections import Counter

# ---------------------------------------------------------------------------------------------------	
def DoOffspringSex(Bearpairs,Femalepercent,CDpairs,equalsexratio):
	'''
	DoOffspringSex()
	This function assigns the sex of each offspring.
	'''	
	# Create empty variable for storing individual offspring information
	offspring=[]
	if equalsexratio == 'WrightFisher':
		# Error statement for even population
		if np.mod(len(Bearpairs),2) != 0:
			print('Special case for Wrigth Fisher assumption, must be even population.')
			sys.exit(-1)
		offsex = np.append(np.zeros(len(Bearpairs)/2,"int"),np.ones(len(Bearpairs)/2,"int"))
		np.random.shuffle(offsex)
		
	# Loop through each mate pair
	for i in range(len(Bearpairs)):
		
		# If equal sex ratio is N
		if equalsexratio == 'N':
		
			# And then loop through each offspring from that mate pair
			for j in range(Bearpairs[i][2]):
				
				# Select sex of the jth offspring - select a random number
				randsex = int(100*np.random.uniform())
				
				# If that random number is less the Femalepercent, assign it to be a female
				if randsex < Femalepercent:
					offsex = 0
				
				# If the random number is greater than the Femalepercent, assign it to be a male
				else:
					offsex = 1

				# And then append all information onto a list storage variable offspring
				offspring.append([Bearpairs[i][0],Bearpairs[i][1],CDpairs[i][0],CDpairs[i][1],offsex])
		
		# If equal sex ratio is Y
		elif equalsexratio == 'AtBirth':
		
			# Error statement for if not even number
			if np.mod(Bearpairs[i][2],2) != 0:
				print('Equal sex ratio does not work for odd litter size. Use offno = 3 and even lambda value.')
				sys.exit(-1)
			
			offsex = np.append(np.zeros(Bearpairs[i][2]/2,"int"),np.ones(Bearpairs[i][2]/2,"int"))
			np.random.shuffle(offsex)
			
			# And then loop through each offspring from that mate pair
			for j in range(Bearpairs[i][2]):
			
				# And then append all information onto a list storage variable offspring
				offspring.append([Bearpairs[i][0],Bearpairs[i][1],CDpairs[i][0],CDpairs[i][1],offsex[j]])
		
		# WRightFisher - assign equal males and females
		elif equalsexratio == 'WrightFisher':
			# And then loop through each offspring from that mate pair
			for j in range(Bearpairs[i][2]):
			
				# And then append all information onto a list storage variable offspring
				offspring.append([Bearpairs[i][0],Bearpairs[i][1],CDpairs[i][0],CDpairs[i][1],offsex[i]])		
	
	# Variables returned
	return offspring
	
	# End::DoOffspringSex()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringRandom(Bearpairs,CDpairs,lmbdavals,age,subpop):
	'''
	DoOffspringRandom()
	This function chooses a random number of 
	offspring for a mated pair.
	'''	
	
	# Loop through each mate pair
	for i in range(len(Bearpairs)):
		
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If females did mate up, then assign random drawn number
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])			
			# Get the subpop female is from - subtract 1 to count from zero
			# Only if more than one Agevars file is given for subpops
			if len(lmbdavals) > 1:
				Fpop = int(subpop[Bearpairs[i][0]])-1
			else:
				Fpop = 0
			if Fage >= len(lmbdavals[Fpop])+1:
				Fage = len(lmbdavals[Fpop])-1
			lmbda = lmbdavals[Fpop][Fage-1]
			
			# Randomly choose a number between 0 and 4
			randkidno = int(round(lmbda)*np.random.uniform())
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(randkidno)
			CDpairs[i].append(randkidno)
	
	# Variables returned
	tupDoOffRand = Bearpairs,CDpairs
	return tupDoOffRand
	# End::DoOffspringRandom()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringPoisson(Bearpairs,CDpairs,lmbdavals,age,subpop):
	'''
	DoOffspringPoisson()
	This function chooses a number of offspring 
	from a Poisson distribution for a mated pair.
	'''		
	
	# Loop through each mate pair
	for i in range(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If female did mate up, then assign a poisson draw with mean lmbda
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])			
			# Get the subpop female is from - subtract 1 to count from zero
			# Only if more than one Agevars file is given for subpops
			if len(lmbdavals) > 1:
				Fpop = int(subpop[Bearpairs[i][0]])-1
			else:
				Fpop = 0
			if Fage >= len(lmbdavals[Fpop])+1:
				Fage = len(lmbdavals[Fpop])-1
			lmbda = lmbdavals[Fpop][Fage-1]
						
			# Poisson number with mean, lambda
			poissonkidno = poisson(lmbda)
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(poissonkidno)
			CDpairs[i].append(poissonkidno)
			
	# Variables returned
	tupDoOffPois = Bearpairs,CDpairs
	return tupDoOffPois
	# End::DoOffspringPoisson()

# ---------------------------------------------------------------------------------------------------	
def DoOffspringConstant(Bearpairs,CDpairs,lmbdavals,age,subpop):
	'''
	DoOffspringConstant()
	This function chooses a constant number of 
	offspring for each mated pair.
	'''	
	
	# Loop through each mate pair
	for i in range(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If females did mate up, then assign lmbda constant
		else:
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			# Get the subpop female is from - subtract 1 to count from zero
			# Only if more than one Agevars file is given for subpops
			if len(lmbdavals) > 1:
				Fpop = int(subpop[Bearpairs[i][0]])-1
			else:
				Fpop = 0
			if Fage >= len(lmbdavals[Fpop])+1:
				Fage = len(lmbdavals[Fpop])-1
			lmbda = lmbdavals[Fpop][Fage-1]
			
			# Assign a 1 [F,M,#offspring]
			Bearpairs[i].append(int(round(lmbda)))
			CDpairs[i].append(int(round(lmbda)))
			
	# Variables returned
	tupDoOffConst = Bearpairs,CDpairs
	return tupDoOffConst
	# End::DoOffspringConstant()
	
# ---------------------------------------------------------------------------------------------------	
def DoOffspringEqual(offspring,lmbdavals,age,subpop):
	'''
	DoOffspringEqual()
	This function chooses lmbda number of offspring for each female.
	'''	
	# import up front not working...temporary.
	from sets import Set
	
	# If offspring list is not empty
	if len(offspring) != 0:		
				
		# Get array and list of females
		tempOff = np.asarray(offspring)
		tempFemales = list(tempOff[:,0])
		
		# Get unique Female locations
		uniqueSet = set(item for item in tempFemales)
		
		# Loop through unique Set of females and create sample list
		returnOff = []
		for item in uniqueSet:
			
			# Get age of unique female
			Fage = int(age[item])
			# Get the subpop female is from - subtract 1 to count from zero
			# Only if more than one Agevars file is given for subpops
			if len(lmbdavals) > 1:
				Fpop = int(subpop[item])-1
			else:
				Fpop = 0
			if Fage >= len(lmbdavals[Fpop])+1:
				Fage = len(lmbdavals[Fpop])-1
			lmbda = lmbdavals[Fpop][Fage-1]
			 
			# Get where this female is in offspring list
			tempFloc = np.where(tempFemales == item)[0]
			# Then random sample lmbda of them
			sampOffloc = np.random.choice(tempFloc,int(round(lmbda)),replace=False).tolist()
			for i in range(len(sampOffloc)):
				returnOff.append(offspring[sampOffloc[i]])
	else:
		returnOff = offspring
		
	return returnOff
	# End::DoOffspringEqual()

# ---------------------------------------------------------------------------------------------------
def DoOffspringNormal(Bearpairs,CDpairs,lmbdavals,sigmavals,age,subpop):
	
	# Loop through each mate pair
	for i in range(len(Bearpairs)):
	
		# If female did not mate up, then assign 0 offspring
		if Bearpairs[i][1] == -9999:
			Bearpairs[i].append(0)
			CDpairs[i].append(0)
		
		# If female did mate up, then assign a normal draw
		else:
			
			# Get age of female
			Fage = int(age[Bearpairs[i][0]])
			# Get the subpop female is from - subtract 1 to count from zero
			# Only if more than one Agevars file is given for subpops
			if len(lmbdavals) > 1:
				Fpop = int(subpop[Bearpairs[i][0]])-1
			else:
				Fpop = 0
			
			if Fage >= len(lmbdavals[Fpop])+1:
				Fage = len(lmbdavals[Fpop])-1
			lmbda = lmbdavals[Fpop][Fage-1]
			sigma = sigmavals[Fpop][Fage-1]
			
			# Normal
			if lmbda != 0:
				if sigma == 0:
					print('Normal draw needs std for offspring number.')
					sys.exit(-1)
				normalkidno = int(round(np.random.normal(lmbda,sigma)))
				if normalkidno < 0:
					normalkidno = 0
			else:
				normalkidno = 0
			
			# Append Offspring number to end of Pairs [F,M,#offspring]
			Bearpairs[i].append(normalkidno)
			CDpairs[i].append(normalkidno)
			
	# Variables returned
	tupDoOff = Bearpairs,CDpairs
	return tupDoOff
	# End::DoOffspringNormal()
	
# ---------------------------------------------------------------------------------------------------	 
def DoCDInfectAndTwinning(offspring,infection,transmissionprob,twinning,Twins):
	'''
	DoCDInfect()
	This function determins whether the offspring gets infected
	or not from the parents infection.  Vertical transmission.
	DoTwinning() 
	Also check for twins and split egg, assinging a unique ID
	'''	
	
	count_twins = 0
	# Get infection status of offspring also looping through 'egg' to see if it splits
	for ioff in range(len(offspring)):
	
		# If parent has infection
		if infection[int(offspring[ioff][0])] == 1 or\
		infection[int(offspring[ioff][1])] == 1:
		
			# Get a random number
			randinfection = np.random.uniform()
			
			# If random flip is less than transmissionprob
			if randinfection < transmissionprob:
			
				# Then append infection status to offspring 
				offspring[ioff].append(1)
				
			# If offspring does not get infected
			else:
				offspring[ioff].append(0)
			
		# If offspring does not get infected.
		else:
		
			# Then append infection status to offspring
			offspring[ioff].append(0)
			
		# Twinning check here - Get a random number and check probability
		randtwin = np.random.uniform()
		if randtwin < twinning:
			# Twinning happens
			offspring[ioff].append('T'+str(ioff)) # Gives unique ID to this twin
			offspring.append(offspring[ioff]) # Then copy this egg				
			count_twins = count_twins + 1	
		else:
			offspring[ioff].append('-9999')
	Twins.append(count_twins)
	
	# Return variables
	return offspring
		
	# End::DoCDInfect()	
	
# ----------------------------------------------------------------------------------------------	 
def DoOffMortality(offspring,Mnewmortperc,Fnewmortperc,equalsexratio,offno,age,lmbda,sex,Track_MOffDeaths,Track_FOffDeaths,subpop):
	'''
	DoOffMortality()
	DoOffMortality()
	Mortality functions for age 0 class.
	'''
	
	if len(offspring) != 0:
				
		# Special case for stable age distribution - adjust survival number
		if offno == '6':
			
			if equalsexratio == 'AtBirth':
				print('Not operating currently with special case of offno=6 for equalsexratio=AtBirth.')
				sys.exit(-1)			
			
			tempoffspring = [] # storage for survived offspring as loop through age classes
			
			if (np.asarray(Mnewmortperc[0]) != np.asarray(Fnewmortperc[0])):
				print('Warning: Age 0 mortality is different for males and females, and option 6 offno is being used. Using average survival.')
			s0 = 1. - ((Mnewmortperc[0] + Fnewmortperc[0]) / 2.)
			s0_pops = 1. - ((np.asarray(Mnewmortperc) + np.asarray(Fnewmortperc)) / 2.)
			
			# Get mothers locations
			motherspots = np.asarray(offspring,dtype = 'int')[:,0]
			motherpops = np.asarray(subpop)[motherspots]
			
			# Loop through each supop
			for ipop in range(len(Mnewmortperc)):			
				
				# Get total in this subpop and their ages for Nt calculation
				thispop = np.where(np.asarray(subpop,dtype='|U10') == str(ipop+1))[0]
				thispop_countage = Counter(np.asarray(age)[thispop]) # string with NAs
				
				# Get the mothers in this subpop, this is their index location back to offspring array
				mothers_inthispop_index = np.where(motherpops == str(ipop+1))[0]
				mothers_inthispop_FID = motherspots[np.where(motherpops == str(ipop+1))[0]]

				# Is there any mothers here?
				if len(mothers_inthispop_index) == 0:
					continue				
				
				# Get the ages of the mothers in this pop
				mothers_ages = np.asarray(np.asarray(age)[mothers_inthispop_FID],dtype = 'int')
				
				# Loop through the age cohorts of mothers and adjust offspring survival
				mothers_countage = Counter(np.asarray(mothers_ages,dtype=np.int8))
				
				# Loop through age classes being careful of indexing
				for iage in range(len(lmbda[0])):
					# Get the age here, 1+
					ageindex = iage + 1				
				
					# Get total individuals in this age class
					Nt = thispop_countage[str(ageindex)]
					# Get total females in this age class
					Ft = mothers_countage[ageindex]			
					# Get number of offspring produced by each age class
					Ft_offno = Ft * lmbda[ipop][iage]
								
					# Females here?
					if Ft == 0:
						continue
					
					# Get corrected survival number
					correctedS0 = s0_pops[ipop]*Nt/Ft
					# Get the number of offspring survivors from this cohort
					nooffsurvivors = int(round(correctedS0))	
					
					# Locate the mothers of this age class, and get original index back
					mothers_thisage = mothers_inthispop_index[np.where(mothers_ages == ageindex)[0]]
					
					# Case where corrected number is more than offspring here
					if nooffsurvivors > len(mothers_thisage): 
						nooffsurvivors = len(mothers_thisage)
						
					# Randomly sample the index for survivors
					offsurvived_index = np.random.choice(mothers_thisage,nooffsurvivors,replace=False).tolist()
					
					# Then store those offspring
					for ioff in range(len(offsurvived_index)):
						tempoffspring.append(offspring[offsurvived_index[ioff]])					
					
			# Tracking
			Track_MOffDeaths.append((len(offspring)-len(tempoffspring))/2.)
			Track_FOffDeaths.append((len(offspring)-len(tempoffspring))/2.)
							
		# Other cases not offno 6
		else:# Get number of survivors for female and male
			if equalsexratio != 'AtBirth':
				tempoffspring = []
				tempTrackMmort = []
				tempTrackFmort = []
				# Skip this loop if mortality is 0, all survive
				if not ((sum(Fnewmortperc)) == 0 and (sum(Mnewmortperc) == 0)):
					# Loop through each offspring
					for ioff in range(len(offspring)):
						thisoff = offspring[ioff]
						# Get subpop from, but only if more than one AgeVars given
						if len(Fnewmortperc) > 1:
							offpop = int(subpop[thisoff[0]])-1 #subract 1 to count from 0
						else:
							offpop = 0
						offsex = thisoff[4]
						if offsex == 0:
							offmort = Fnewmortperc[offpop]
						else:
							offmort = Mnewmortperc[offpop]
						# See if survives
						randno = np.random.uniform()
						if randno < offmort: # mortality occurs
							# Tracking
							if offsex == 0:
								tempTrackFmort.append(1)
							else:
								tempTrackMmort.append(1)
						else: # offspring survives
							tempoffspring.append(thisoff)
					
					# Tracking
					Track_MOffDeaths.append(sum(tempTrackMmort))
					Track_FOffDeaths.append(sum(tempTrackFmort))
				else:
					# Tracking
					Track_MOffDeaths.append(0)
					Track_FOffDeaths.append(0)
					tempoffspring = offspring
			# For AtBrith option 
			else:
				if len(Mnewmortperc) > 1:
					print('Multiple AgeVars files (e.g., subpops) not currently operating with AtBirth equal sex ratio option.')
					sys.exit(-1)
				else:
					offsex = np.asarray(offspring)[:,4]
					Foffsex = np.where(offsex == '0')[0]
					Moffsex = np.where(offsex == '1')[0]
					# Index location
					Foffdeaths = int(round((Fnewmortperc[0])*len(Foffsex)))
					Moffdeaths = int(round((Mnewmortperc[0])*len(Moffsex)))			
					
					# Shuffle Female and male index list
					shuffle(Foffsex)
					shuffle(Moffsex)
					
					if (Foffdeaths + Moffdeaths) > 0:
						print('Warning: Equal sex ratio AtBirth is specified with different age 0 mortality values. Using total age 0 deaths.')
					offdeaths = Foffdeaths + Moffdeaths
					
					# Remove equally from F and M 
					tempFoff = Foffsex[0:len(Foffsex)-(offdeaths/2)]
					tempMoff = Moffsex[0:len(Moffsex)-(offdeaths/2)]					
					
					# Get index death locations for F and M
					delFoff = Foffsex[0:offdeaths/2]
					delMoff = Moffsex[0:offdeaths/2]
					deloff = np.concatenate((delFoff,delMoff),axis=0)
					
					for ioff in sorted(deloff,reverse=True):
						del offspring[ioff]
					tempoffspring = offspring
					# Tracking
					Track_MOffDeaths.append(Moffdeaths)
					Track_FOffDeaths.append(Foffdeaths)
	
	# Else no offspring
	else:
		Track_MOffDeaths.append(0)
		Track_FOffDeaths.append(0)
		tempoffspring = offspring
		
	return tempoffspring
	
	# End::DoOffMortality()
	
# ---------------------------------------------------------------------------------------------------	 
def DoOffspring(offno,lmbda,Bearpairs,CDpairs,\
Femalepercent,Births,infection,transmissionprob,equalsexratio,Mnewmortperc,Fnewmortperc,Track_MOffDeaths,Track_FOffDeaths,sigma,age,sex,twinning,Twins,subpop):
	'''
	DoOffspring()
	Choose numberof Offspring for each mated pair
	Input: selection choice for offspring number distrubution draw
	offno, Bearpairs, lmbda.
	Output: Bear Pairs + # of offspring [Female,Male,#Offspring]
	Also add # offspring to CDpairs
	offspring returned [Female,Male,Female,Male,sex,infection,TwinningID]
	'''
	
	# Function 1 is a uniform random draw between 0 and lmdba number	
	if (offno=='1'):
		
		tupDoOffRand = DoOffspringRandom(Bearpairs,CDpairs,lmbda,age,subpop)
		
		Bearpairs = tupDoOffRand[0]
		CDpairs = tupDoOffRand[1]
	
	# Function 2 is a Poisson draw
	elif (offno=='2'):
	
		tupDoOffPois = DoOffspringPoisson(Bearpairs,CDpairs,lmbda,age,subpop)
		
		Bearpairs = tupDoOffPois[0]
		CDpairs = tupDoOffPois[1]
		
	# Function 3 and 4 is a constant of lmbda offspring per each pairing
	elif (offno=='3' or offno=='4' or offno == '6'):
	
		tupDoOffConst = DoOffspringConstant(Bearpairs,CDpairs,\
		lmbda,age,subpop)
		
		Bearpairs = tupDoOffConst[0]
		CDpairs = tupDoOffConst[1]
	
	# Function 5 is a normal draw
	elif (offno=='5'):
		
		tupDoOff = DoOffspringNormal(Bearpairs,CDpairs,lmbda,sigma,age,subpop)
		Bearpairs = tupDoOff[0]
		CDpairs = tupDoOff[1]
	
	# Other functions
	else:
		print('This offspring birth rate option (offno) does not exist.')
		sys.exit(-1)
	
	# Assign sex of each offspring
	offspring = DoOffspringSex(Bearpairs,Femalepercent,CDpairs,equalsexratio)
	
	# If function 4 is choosen, then weed the offspring here
	if (offno == '4'):
		
		offspring = DoOffspringEqual(offspring,lmbda,age,subpop)	
	
	# Assign infection to each offspring and Apply twinning possibility here
	offspring = DoCDInfectAndTwinning(offspring,infection,transmissionprob,twinning,Twins)
	
	# Store number of Births
	Births.append(len(offspring))	
	
	# Apply Mortality to the egg class here
	offspring = DoOffMortality(offspring,Mnewmortperc,Fnewmortperc,equalsexratio,offno,age,lmbda,sex,Track_MOffDeaths,Track_FOffDeaths,subpop)
	
	# Return variables from this argument, offspring values are all string now.
	tupDoOff = offspring,len(offspring)
	return tupDoOff
		
	# End::DoOffspring()