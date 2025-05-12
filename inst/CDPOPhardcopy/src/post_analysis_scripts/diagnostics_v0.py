# -----------------------------------------------------------------------------
# diagnostics.py
# 2014 01 29: Erin Landguth
# This script grabs output.csv files, splits bars, and grabs diagnostics for checking CDPOP/CDFISH/CDmetaPOP
# v0 - Initial script
# ----------------------------------------------------------------------------- 

# Load modules
import os,pdb,pandas
from pylab import *	
import scipy as sp			

# Numpy functions
try:
	import numpy as np 
	from numpy.random import *
except ImportError:
	raise ImportError("Numpy required.")

# ---------
# User info
# ---------
# Directory locations of folders
#dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanStable_expN_120gens_2K500N_Kenv1000_1411750885/'
#plottitle = 'PAN Stable 2K Exponential N'
dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanStable_expAtBirth_120gens_2K500N_Kenv1000_1411751507/'
plottitle = 'PAN Stable 2K Exponential AtBirth'
dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanStable_logAtBirth_120gens_2K500N_Kenv1000_1411751520/'
plottitle = 'PAN Stable 2K Logistic AtBirth'
#dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanStable_logN_120gens_2K500N_Kenv1000_1411748653/'
#plottitle = 'PAN Stable 2K Logistic N'
dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanIncrease_exp_50gens_2K500N_1411754599/'
plottitle = 'PAN Increasing 2K Exponential N'
dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanStable_exps0pt26_120gens_2K500N_1411753815/'
plottitle = 'PAN Increasing s0=0.26 2K Exponential N'

dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanDecrease_exp_100gens_2K500N_1411759172/'
plottitle = 'PAN Decreasing s0=0.24 2K Exponential N'

dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/IBDStable_exp_400gens_2K500N_testHe_1411763780/'
plottitle = 'IBDMax Stable s0=0.25... 2K Exponential N'
dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanFluctuate_exp_94gens_2K500N_1411764557/'
plottitle = 'PAN Fluctuate s0=0.25... 2K Exponential N'

#dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/IBDIncrease_pt29_exp_50gens_2K500N_1411766939/'
#plottitle = 'IBD Increase s0=0.29 2K Exponential N'

#dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/IBDDecrease_pt21_exp_50gens_2K500N_1411766969/'
#plottitle = 'IBD Decrease s0~0.21 2K Exponential N'

#dir = 'D:/projects/CDPOP/TemporalDynamics/Runs/2Kdata/PanFluctuate_log_94gens_2K500N_1411831887/'
#plottitle = 'PAN Fluctuate s0~0.29 to 0.21 2K Logistic N'

savename = ''
outdir = dir

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 94 # Number of years 
nthfile = list(range(0,gen,1))
maxA = 900 # loci * alleles
batchno = 1 # Number of batches
mcno = 7 # Number of MCs
plottime = np.asarray([1,9,19])
K = 2000

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# ---------------------------------
# Storage variables
# ---------------------------------
# Store values
N = [] # 1
N_age = []
growthrate = []
totF = [] # 4
totM = [] #5
breedF = [] #6
breedM = [] #8
births = [] #13
mort_age = [] # 15

alleles = [] #16
He = [] #17
Ho = [] #18
'''
outputtitle = ['Year','Population','Population_Age','ToTFemales','ToTMales',\
	'BreedFemales','BreedFemales_Age','BreedMales','BreedEvents_Females','Females_NoMate',\
	'Migrants','DisperseDeaths','Births','EggDeaths',\
	'AgeDeaths','Alleles','He','Ho','Mutations',\
	'MateDistED','MateDistEDstd','Female_DispDistED','Female_DispDistEDstd',\
	'Male_DispDistED','Male_DispDistEDstd','MateDistCD','MateDistCDstd',\
	'Female_DispDistCD','Female_DispDistCDstd','Male_DispDistCD','Male_DispDistCDstd',\
	'p1','p2','q1','q2','Infected','SubpopImmigration','SubpopEmigration','FemalesMeanMate',\
	'MalesMeanMate','FemalesSDMate','MalesSDMate','OpenLocations','CouldNotDisperse']
'''
# Loop through batches
for ibatch in range(batchno):

	# Add storage spot
	N.append([])
	N_age.append([])
	growthrate.append([])
	totM.append([])
	totF.append([])
	breedF.append([])
	breedM.append([])
	births.append([])	
	mort_age.append([]) #11 with total
	alleles.append([])
	He.append([])
	Ho.append([])
	
	# Loop through MCs
	for imc in range(mcno):
	
		# Open output.csv file in folder
		inputfile = open(dir+'batchrun'+str(ibatch)+'mcrun'+str(imc)+'/output.csv')
		
		# Read lines from the file
		lines = inputfile.readlines()

		#Close the file
		inputfile.close()

		# Create an empty matrix to append to
		values = []

		# Split up each line in file and append to empty matrix for generation specified
		for i in range(len(lines)):
			thisline = lines[i].split(',')
			values.append(thisline)

		# Delete lines
		del(lines)
		
		# Store values
		N[ibatch].append([])
		N_age[ibatch].append([])
		growthrate[ibatch].append([])
		totM[ibatch].append([])
		totF[ibatch].append([])
		breedF[ibatch].append([])
		breedM[ibatch].append([])
		births[ibatch].append([])
		mort_age[ibatch].append([]) #11 with total
		alleles[ibatch].append([])
		He[ibatch].append([])
		Ho[ibatch].append([])
		
		# Then Loop through generations/time
		for iout in range(gen):
		
			alleles[ibatch][imc].append([])
			He[ibatch][imc].append([])
			Ho[ibatch][imc].append([])
			mort_age[ibatch][imc].append([])
			N_age[ibatch][imc].append([])
			
			# Patch split - grab the first value
			N[ibatch][imc].append(float(int(values[1+iout][1].split('|')[0])))
			totF[ibatch][imc].append(int(values[1+iout][4].split('|')[0]))
			totM[ibatch][imc].append(int(values[1+iout][5].split('|')[0]))
			breedF[ibatch][imc].append(int(values[1+iout][6].split('|')[0]))
			breedM[ibatch][imc].append(int(values[1+iout][8].split('|')[0]))
			births[ibatch][imc].append(int(values[1+iout][13].split('|')[0]))
			if values[1+iout][3] == 'NA':
				growthrate[ibatch][imc].append(np.nan)
			else:
				growthrate[ibatch][imc].append(float(values[1+iout][3]))
						
			# Age split
			for j in range(len(values[1+iout][15].split('|'))-1):
				mort_age[ibatch][imc][iout].append(float(values[1+iout][15].split('|')[j]))
				N_age[ibatch][imc][iout].append(float(values[1+iout][2].split('|')[j]))				
						
			# Grab all patch values - patch values with total
			for j in range(len(values[1+iout][16].split('|'))-1):
				alleles[ibatch][imc][iout].append(float(values[1+iout][16].split('|')[j]))
				He[ibatch][imc][iout].append(float(values[1+iout][17].split('|')[j]))
				Ho[ibatch][imc][iout].append(float(values[1+iout][18].split('|')[j]))
							
			
# Turn into arrays
births = np.asarray(births)
totF = np.asarray(totF)
alleles = np.asarray(alleles)
He = np.asarray(He)
Ho = np.asarray(Ho)
allD = (alleles - 1) / (maxA - 1.)
N = np.asarray(N)
totM = np.asarray(totM)
breedF = np.asarray(breedF)
breedM = np.asarray(breedM)
mort_age = np.asarray(mort_age)
N_age = np.asarray(N_age)
growthrate = np.asarray(growthrate)

# --------------------------------------------
# Get mean over Monte Carlosfor each batch run
# --------------------------------------------
# Get Mean, SD, error, Left and Right
N_m = np.nansum(N[:][:],axis=1)/mcno
N_sd = np.std(N[:][:],axis=1)	
error = qnorm*N_sd/(mcno)
N_l = N_m-error
N_r = 	N_m+error
N_min = np.min(N[:][:],axis=1)
N_max = np.max(N[:][:],axis=1)	

growthrate_m = np.nansum(growthrate[:][:],axis=1)/mcno
growthrate_sd = np.std(growthrate[:][:],axis=1)	
error = qnorm*growthrate_sd/(mcno)
growthrate_l = growthrate_m-error
growthrate_r = 	growthrate_m+error
growthrate_min = np.min(growthrate[:][:],axis=1)
growthrate_max = np.max(growthrate[:][:],axis=1)

totM_m = np.nansum(totM[:][:],axis=1)/mcno
totM_sd = np.std(totM[:][:],axis=1)	
error = totM_sd*qnorm/mcno
totM_l = totM_m-error
totM_r = 	totM_m+error
totM_min = np.min(totM[:][:],axis=1)
totM_max = np.max(totM[:][:],axis=1)

totF_m = np.nansum(totF[:][:],axis=1)/mcno
totF_sd = np.std(totF[:][:],axis=1)	
error = qnorm*totF_sd/(mcno)
totF_l = totF_m-error
totF_r = 	totF_m+error
totF_min = np.min(totF[:][:],axis=1)
totF_max = np.max(totF[:][:],axis=1)

breedF_m = np.nansum(breedF[:][:],axis=1)/mcno
breedF_sd = np.std(breedF[:][:],axis=1)	
error = qnorm*breedF_sd/(mcno)
breedF_l = breedF_m-error
breedF_r = 	breedF_m+error
breedF_min = np.min(breedF[:][:],axis=1)
breedF_max = np.max(breedF[:][:],axis=1)

breedM_m = np.nansum(breedM[:][:],axis=1)/mcno
breedM_sd = np.std(breedM[:][:],axis=1)	
error = qnorm*breedM_sd/(mcno)
breedM_l = breedM_m-error
breedM_r = 	breedM_m+error
breedM_min = np.min(breedM[:][:],axis=1)
breedM_max = np.max(breedM[:][:],axis=1)
births_m = np.nansum(births[:][:],axis=1)/mcno
births_sd = np.std(births[:][:],axis=1)	
error = qnorm*births_sd/(mcno)
births_l = births_m-error
births_r = 	births_m+error
births_min = np.min(births[:][:],axis=1)
births_max = np.max(births[:][:],axis=1)

mort_age_m = np.nansum(mort_age[:][:],axis=1)/mcno
mort_age_sd = np.std(mort_age[:][:],axis=1)	
error = qnorm*mort_age_sd/(mcno)
mort_age_l = mort_age_m-error
mort_age_r = 	mort_age_m+error
mort_age_min = np.min(mort_age[:][:],axis=1)
mort_age_max = np.max(mort_age[:][:],axis=1)

N_age_m = np.nansum(N_age[:][:],axis=1)/mcno
N_age_sd = np.std(N_age[:][:],axis=1)	
error = qnorm*N_age_sd/(mcno)
N_age_l = N_age_m-error
N_age_r = 	N_age_m+error
N_age_min = np.min(N_age[:][:],axis=1)
N_age_max = np.max(N_age[:][:],axis=1)

He_m = np.nansum(He[:][:],axis=1)/mcno
He_sd = np.std(He[:][:],axis=1)	
error = qnorm*He_sd/(mcno)
He_l = He_m-error
He_r = 	He_m+error
He_min = np.min(He[:][:],axis=1)
He_max = np.max(He[:][:],axis=1)

Ho_m = np.nansum(Ho[:][:],axis=1)/mcno
Ho_sd = np.std(Ho[:][:],axis=1)	
error = qnorm*Ho_sd/(mcno)
Ho_l = Ho_m-error
Ho_r = 	Ho_m+error
Ho_min = np.min(Ho[:][:],axis=1)
Ho_max = np.max(Ho[:][:],axis=1)

allD_m = np.nansum(allD[:][:],axis=1)/mcno
allD_sd = np.std(allD[:][:],axis=1)	
error = qnorm*allD_sd/(mcno)
allD_l = allD_m-error
allD_r = 	allD_m+error

# --------------------------------
# Other summary data 
# -------------------------------- 

# Growth rate - total population Nt / Nt-1
Nt = N[:,:,1:gen] / N[:,:,0:gen-1]
N_growth_pop_m = np.nansum(Nt[:][:],axis=1)/(mcno)
N_growth_pop_sd = np.nanstd(Nt[:][:],axis=1)	
error = qnorm*N_growth_pop_sd/(mcno)
N_growth_pop_l = N_growth_pop_m-error
N_growth_pop_r = 	N_growth_pop_m+error
N_growth_pop_min = np.nanmin(Nt[:][:],axis=2)
N_growth_pop_max = np.nanmax(Nt[:][:],axis=2)

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot all initial N - total population - fill between
# ----------------------------------------------------
figure()
plot(nthfile,N_m[0],'k-')
plot(nthfile,N_l[0],'r.-', ms = 10)
plot(nthfile,N_r[0],'r.-', ms = 10)
fill_between(nthfile, N_min[0], N_max[0],color = 'b')
fill_between(nthfile, N_l[0], N_r[0],color = 'r')
xlabel('Time',fontsize=18)
ylabel('Population Total',fontsize=18)
title(plottitle+' R='+str(round(np.nanmean(growthrate_m[0]),4)),fontsize=16)
axis([-0.1,gen,0,K])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+savename+'N.png',dpi=savedpi)


# Plot All N as bar - stacked side by side
# ------------------------------------------
ind = np.asarray(nthfile)  # the x locations for the groups
width = 0.15       # the width of the bars
fig, ax = plt.subplots()

rects1 = ax.bar(ind[0:len(plottime)], N_m[0][plottime], width, color='grey', edgecolor='black', hatch="/")
yminus = N_m[0][plottime]-N_min[0][plottime]
yplus =  N_m[0][plottime]-N_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+width/2., N_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects2 = ax.bar(ind[0:len(plottime)]+width, totF_m[0][plottime], width, color='grey', edgecolor='black', hatch="//")
yminus = totF_m[0][plottime]-totF_min[0][plottime]
yplus =  totF_m[0][plottime]-totF_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+3*width/2., totF_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects3 = ax.bar(ind[0:len(plottime)]+2*width, totM_m[0][plottime], width, color='grey', edgecolor='black')
yminus = totM_m[0][plottime]-totM_min[0][plottime]
yplus =  totM_m[0][plottime]-totM_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+5*width/2., totM_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects4 = ax.bar(ind[0:len(plottime)]+3*width, breedF_m[0][plottime], width, color='grey', edgecolor='black', hatch="x")
yminus = breedF_m[0][plottime]-breedF_min[0][plottime]
yplus =  breedF_m[0][plottime]-breedF_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+7*width/2., breedF_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects5 = ax.bar(ind[0:len(plottime)]+4*width, breedM_m[0][plottime], width, color='grey', edgecolor='black', hatch="\\")
yminus = breedM_m[0][plottime]-breedM_min[0][plottime]
yplus =  breedM_m[0][plottime]-breedM_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+9*width/2., breedM_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects6 = ax.bar(ind[0:len(plottime)]+5*width, births_m[0][plottime], width, color='grey', edgecolor='black', hatch="+")
yminus = births_m[0][plottime]-births_min[0][plottime]
yplus =  births_m[0][plottime]-births_max[0][plottime]
ax.errorbar(ind[0:len(plottime)]+11*width/2., births_m[0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

# add some
ax.set_ylabel('N values',fontsize=18)
ax.set_xlabel('Time',fontsize=18)
ax.set_title(plottitle,fontsize=21)
ax.set_xticks(ind[0:len(plottime)]+2*width)
ax.set_xticklabels( (np.asarray(plottime,dtype='str')),ha='center' )
ax.legend( (rects1[0], rects2[0],rects3[0],rects4[0],rects5[0],rects6[0]), ('N', 'F','M','Breed F','Breed M','Births') )
ax.axis([-0,len(plottime),0,K])

# Plot age N
# ----------

colors = ['k','b','r','y','g','m']
tick = ['k-','b.-','r--','y,','g,-','m-']
leg = ['0','1','2','3','4','5']
figure()
plot(nthfile,births_m[0],'c-',label='eggs')
fill_between(nthfile, births_l[0], births_r[0],color = 'c')
for i in range(len(N_age_m[0][0])):
	plot(nthfile,N_age_m[0][:,i],tick[i],label=leg[i])
	fill_between(nthfile, N_age_l[0][:,i], N_age_r[0][:,i],color=colors[i])
#plot(nthfile,N_l[0],'r.-', ms = 10)
#plot(nthfile,N_r[0],'r.-', ms = 10)
xlabel('Time',fontsize=18)
ylabel('Population Age Total',fontsize=18)
title(plottitle,fontsize=21)
axis([-0.1,gen,0,K])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+savename+'N_age.png',dpi=savedpi)

# Plot He Ho
# ----------------------------------------------------
figure()
plot(nthfile,He_m[0][:,0],'k-',label='He')
plot(nthfile,Ho_m[0][:,0],'r.',label='Ho')
#fill_between(nthfile, He_l[0][:,0], He_r[0][:,0],color = 'b')
#fill_between(nthfile, Ho_l[0][:,0], Ho_r[0][:,0],color = 'r')
xlabel('Time',fontsize=18)
ylabel('Heterozygosity',fontsize=18)
title(plottitle+' R='+str(round(np.mean(growthrate_m[0]),4)),fontsize=16)
axis([-0.1,gen,0,1.0])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
savefig(dir+savename+'H.png',dpi=savedpi)


'''
# Plot all mortality - stacked side by side
# ------------------------------------------
ind = np.asarray(nthfile)  # the x locations for the groups
width = 0.13       # the width of the bars
fig, ax = plt.subplots()

rects1 = ax.bar(ind[0:len(plottime)], mort_age_m[0][:,0][plottime], width, color='grey', edgecolor='black', hatch="/")
yminus = mort_age_m[0][:,0][plottime]-mort_age_min[0][:,0][plottime]
yplus =  mort_age_m[0][:,0][plottime]-mort_age_max[0][:,0][plottime]
ax.errorbar(ind[0:len(plottime)]+width/2., mort_age_m[0][:,0][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

rects2 = ax.bar(ind[0:len(plottime)]+width, mort_age_m[0][:,1][plottime], width, color='grey', edgecolor='black', hatch="//")
yminus = mort_age_m[0][:,1][plottime]-mort_age_m[0][:,1][plottime]
yplus =  mort_age_m[0][:,1][plottime]-mort_age_m[0][:,1][plottime]
ax.errorbar(ind[0:len(plottime)]+3*width/2., mort_age_m[0][:,1][plottime], yerr=[yminus, yplus],fmt='bo',capthick=2)

# add some
ax.set_ylabel('Mortality',fontsize=18)
ax.set_xlabel('Time',fontsize=18)
ax.set_title(plottitle,fontsize=21)
ax.set_xticks(ind[0:len(plottime)]+2*width)
ax.set_xticklabels( (np.asarray(plottime,dtype='str')),ha='center' )
ax.legend( (rects1[0], rects2[0]), ('Age0', 'Age1') )
ax.axis([-0,len(plottime),0,K])
savefig(dir+savename+'Allmortality.png',dpi=savedpi)
'''

show()

