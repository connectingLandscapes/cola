# -----------------------------------------------------------------------------
# plotN.py
# 2017 July: Erin Landguth
# This script grabs output.csv files, splits bars, and plot population total
# and by subpopulation
# v0 - Initial script
# v1 - Option to plot each subpop individually
# ----------------------------------------------------------------------------- 

# Load modules
import os,pdb
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
dir = 'D:/projects/CDPOP/BighornSheep_KimAndrewsIGFP/cdpop_runs/data_8pops/Scenario3_ChangingMortRates_1501692857/'
plottitle = ''

savename = 'ChangingMortRates_'
outdir = dir

savedpi = 300
qnorm = 1.959964 # For CIs, not in function 
gen = 51 # Number of years/gens total 
batchno = 1 # Total Number of batches
mcno = 10 # Number of MCs
plottime = [50] # For plotting age distribution for total population
linemarks = ['-','--','-.','-^','-s','-1','-2','-*','-+']
RULabels = ['37-50','18','58','51','11','30A','37A','13']
ploteachSubPop = True

# -------------
# Initial Setup
# -------------
nthfile = list(range(0,gen,1))
nthfile_labels = list(range(1969,gen+1969,1))

# List folders in this directory
def listdirs(folder):
    return [d for d in (os.path.join(folder, d1) for d1 in os.listdir(folder)) if os.path.isdir(d)]
folderList = listdirs(dir)

# Storage variables
N = [] # 1
N_age = [] #2
growthrate = []#3
totF = [] # 4
totM = [] #5
breedF = [] #6
breedM = [] #8

# -------------------
# Begin batch loop
# -------------------
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
				
		# Then Loop through generations/time
		for iout in range(gen):
		
			# Store values
			N[ibatch][imc].append([])
			N_age[ibatch][imc].append([])
			growthrate[ibatch][imc].append([])
			totM[ibatch][imc].append([])
			totF[ibatch][imc].append([])
			breedF[ibatch][imc].append([])
			breedM[ibatch][imc].append([])
			
			# Subpopulation split - remove total 
			for j in range(1,len(values[1+iout][1].split('|'))-1):			
				N[ibatch][imc][iout].append(float(int(values[1+iout][1].split('|')[j])))
				totF[ibatch][imc][iout].append(int(values[1+iout][4].split('|')[j]))
				totM[ibatch][imc][iout].append(int(values[1+iout][5].split('|')[j]))
				breedF[ibatch][imc][iout].append(int(values[1+iout][6].split('|')[j]))
				breedM[ibatch][imc][iout].append(int(values[1+iout][8].split('|')[j]))
			
			# Single values
			if values[1+iout][3] == 'NA':
				growthrate[ibatch][imc][iout].append(np.nan)
			else:
				growthrate[ibatch][imc][iout].append(float(values[1+iout][3]))
						
			# Age split
			for j in range(len(values[1+iout][2].split('|'))-1):
				N_age[ibatch][imc][iout].append(float(values[1+iout][2].split('|')[j]))				
										
			
# Turn into arrays
totF = np.asarray(totF)
N = np.asarray(N)
totM = np.asarray(totM)
breedF = np.asarray(breedF)
breedM = np.asarray(breedM)
N_age = np.asarray(N_age)
growthrate = np.asarray(growthrate)

# ---------------------------------------------------
# Get summaries for each - total, mean, sd, max, min
# --------------------------------------------------
# Get Mean, SD, error, Left and Right
N_m = np.nansum(N,axis=1)/mcno
N_sd = np.std(N,axis=1)	
error = qnorm*N_sd/(mcno)
N_l = N_m-error
N_r = 	N_m+error
N_min = np.min(N,axis=1)
N_max = np.max(N,axis=1)
# Total
N_total = np.nansum(N,axis=3)
N_total_m = np.nansum(N_total,axis=1)/mcno
N_total_sd = np.std(N_total,axis=1)
error = qnorm*N_total_sd/(mcno)
N_total_l = N_total_m-error
N_total_r = 	N_total_m+error

growthrate_m = np.nansum(growthrate,axis=1)/mcno
growthrate_sd = np.std(growthrate,axis=1)	
error = qnorm*growthrate_sd/(mcno)
growthrate_l = growthrate_m-error
growthrate_r = 	growthrate_m+error
growthrate_min = np.min(growthrate,axis=1)
growthrate_max = np.max(growthrate,axis=1)

totM_m = np.nansum(totM,axis=1)/mcno
totM_sd = np.std(totM,axis=1)	
error = totM_sd*qnorm/mcno
totM_l = totM_m-error
totM_r = 	totM_m+error
totM_min = np.min(totM,axis=1)
totM_max = np.max(totM,axis=1)
# Total
totM_total = np.nansum(totM,axis=3)
totM_total_m = np.nansum(totM_total,axis=1)/mcno
totM_total_sd = np.std(totM_total,axis=1)
error = qnorm*totM_total_sd/(mcno)
totM_total_l = totM_total_m-error
totM_total_r = 	totM_total_m+error

totF_m = np.nansum(totF,axis=1)/mcno
totF_sd = np.std(totF,axis=1)	
error = qnorm*totF_sd/(mcno)
totF_l = totF_m-error
totF_r = 	totF_m+error
totF_min = np.min(totF,axis=1)
totF_max = np.max(totF,axis=1)
# Total
totF_total = np.nansum(totF,axis=3)
totF_total_m = np.nansum(totF_total,axis=1)/mcno
totF_total_sd = np.std(totF_total,axis=1)
error = qnorm*totF_total_sd/(mcno)
totF_total_l = totF_total_m-error
totF_total_r = 	totF_total_m+error

breedF_m = np.nansum(breedF,axis=1)/mcno
breedF_sd = np.std(breedF,axis=1)	
error = qnorm*breedF_sd/(mcno)
breedF_l = breedF_m-error
breedF_r = 	breedF_m+error
breedF_min = np.min(breedF,axis=1)
breedF_max = np.max(breedF,axis=1)
# Total
breedF_total = np.nansum(breedF,axis=3)
breedF_total_m = np.nansum(breedF_total,axis=1)/mcno
breedF_total_sd = np.std(breedF_total,axis=1)
error = qnorm*breedF_total_sd/(mcno)
breedF_total_l = breedF_total_m-error
breedF_total_r = 	breedF_total_m+error

breedM_m = np.nansum(breedM,axis=1)/mcno
breedM_sd = np.std(breedM,axis=1)	
error = qnorm*breedM_sd/(mcno)
breedM_l = breedM_m-error
breedM_r = 	breedM_m+error
breedM_min = np.min(breedM,axis=1)
breedM_max = np.max(breedM,axis=1)
# Total
breedM_total = np.nansum(breedM,axis=3)
breedM_total_m = np.nansum(breedM_total,axis=1)/mcno
breedM_total_sd = np.std(breedM_total,axis=1)
error = qnorm*breedM_total_sd/(mcno)
breedM_total_l = breedM_total_m-error
breedM_total_r = 	breedM_total_m+error

N_age_m = np.nansum(N_age,axis=1)/mcno
N_age_sd = np.std(N_age,axis=1)	
error = qnorm*N_age_sd/(mcno)
N_age_l = N_age_m-error
N_age_r = 	N_age_m+error
N_age_min = np.min(N_age,axis=1)
N_age_max = np.max(N_age,axis=1)

# --------------------------------
# Other summary data 
# -------------------------------- 

# Growth rate - total population Nt / Nt-1
Nt = N[:,:,1:gen] / N[:,:,0:gen-1]
N_growth_pop_m = np.nansum(Nt,axis=1)/(mcno)
N_growth_pop_sd = np.nanstd(Nt,axis=1)	
error = qnorm*N_growth_pop_sd/(mcno)
N_growth_pop_l = N_growth_pop_m-error
N_growth_pop_r = 	N_growth_pop_m+error
N_growth_pop_min = np.nanmin(Nt,axis=2)
N_growth_pop_max = np.nanmax(Nt,axis=2)

# --------------------------------------------------------
# Plotting
# --------------------------------------------------------

# Plot all N - total population - fill between
# ----------------------------------------------------
figure(1)
plot(nthfile,N_total_m[0],'k-')
plot(nthfile,N_total_l[0],'r.-', ms = 10)
plot(nthfile,N_total_r[0],'r.-', ms = 10)
#fill_between(nthfile, N_total_min[0], N_total_max[0],color = 'b')
fill_between(nthfile, N_total_l[0], N_total_r[0],color = 'r')
xlabel('Time',fontsize=18)
ylabel('Population Total',fontsize=18)
title(plottitle+' R='+str(round(np.nanmean(growthrate_m[0]),4)),fontsize=16)
#axis([-0.1,gen,0,K])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
xticks(nthfile,nthfile_labels,rotation='vertical')
#savefig(dir+savename+'N_total.png',dpi=savedpi)

# Plot all subpop Ns
# ------------------
figure(2)
for ipop in range(len(N_m[0][0])):
	plot(nthfile,N_m[0,:,ipop],linemarks[ipop],label=RULabels[ipop])
xlabel('Time',fontsize=18)
ylabel('RU Totals',fontsize=18)
#title(plottitle+' R='+str(round(np.nanmean(growthrate_m[0]),4)),fontsize=16)
#axis([-0.1,gen,0,K])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
xticks(nthfile,nthfile_labels,rotation='vertical')
#savefig(dir+savename+'RU_total.png',dpi=savedpi)

# Plot Total Females and Males
# --------------------------------------------- 
figure(3)
plot(nthfile,totM_total_m[0],'b-',label='Total Males')
plot(nthfile,totF_total_m[0],'r-',label='Total Females')
plot(nthfile,breedM_total_m[0],'b-.',label='Breed Males')
plot(nthfile,breedF_total_m[0],'r-.',label='Breed Females')
#plot(nthfile,totM_total_l[0],'r.-', ms = 10)
#plot(nthfile,totM_total_r[0],'r.-', ms = 10)
#fill_between(nthfile, N_total_min[0], N_total_max[0],color = 'b')
#fill_between(nthfile, N_total_l[0], N_total_r[0],color = 'r')
xlabel('Time',fontsize=18)
ylabel('Sex Totals',fontsize=18)
#title(plottitle+' R='+str(round(np.nanmean(growthrate_m[0]),4)),fontsize=16)
#axis([-0.1,gen,0,K])
legend(loc=0)
# Updating fontsize on axes
fontsize=19
ax = gca()
for tick in ax.xaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
	tick.label1.set_fontsize(fontsize)
xticks(nthfile,nthfile_labels,rotation='vertical')
#savefig(dir+savename+'Sex_totals.png',dpi=savedpi)

# Plot N initial - age class
# -----------------------------
ind = np.arange(len(N_age_m[0][0])) # the x locations for the groups
width = 0.15                      # the width of the bars

# Loop through each year to plot
for it in plottime:
	
	fig = plt.figure(4)
	ax = fig.add_subplot(111)
	
	rects0 = ax.bar(ind,N_age_m[0][it],width,color='red',edgecolor='black',hatch="/",yerr=N_age_sd[0][it],error_kw=dict(elinewidth=2,ecolor='black'))		
		
	# axes and labels
	ax.set_ylabel('N',fontsize=18)
	ax.set_xlabel('Age',fontsize=18)
	ax.set_title(plottitle+'Age N at Year ' +str(nthfile_labels[it]),fontsize=21)
	ax.set_xlim(-width,len(ind)+width)
	#ax.set_ylim(0, max(np.max(N_init_age_m[0][it]),np.max(N_init_age_m[1][it]),np.max(N_init_age_m[2][it]),np.max(N_init_age_m[3][it]),np.max(N_init_age_m[4][it])))
	xTickMarks = [str(i) for i in range(len(N_age_m[0][0]))]
	ax.set_xticks(ind+width)
	xtickNames = ax.set_xticklabels(xTickMarks)
	plt.setp(xtickNames, rotation=0)
	#ax.legend((rects0[0]),'Year '+str(nthfile_labels[it]),loc=0)
	savefig(dir+savename+'AGEDistribution_Year'+str(it)+'.png',dpi=savedpi)

# Plot each subpop separately
# ---------------------------
if ploteachSubPop:
	for ipop in range(len(N_m[0][0])):
		figure(5+ipop)
		plot(nthfile,N_m[0,:,ipop],linemarks[ipop])
		xlabel('Time',fontsize=18)
		ylabel('Total N',fontsize=18)
		title('RU '+RULabels[ipop],fontsize=16)
		# Updating fontsize on axes
		fontsize=19
		ax = gca()
		for tick in ax.xaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
		for tick in ax.yaxis.get_major_ticks():
			tick.label1.set_fontsize(fontsize)
		xticks(nthfile,nthfile_labels,rotation='vertical')
		savefig(dir+savename+'RU_'+RULabels[ipop]+'.png',dpi=savedpi)
	
show()
