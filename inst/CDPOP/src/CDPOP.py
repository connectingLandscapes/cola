# CDPOP.py
# Author: Erin L Landguth
# Created: February 2008
# v 1.2 Release: December 2011
# ----------------------------------------------------------------------------
# General CDPOP information
appName = "CDPOP"
appVers = "version 1.3.20"
appRele = "2024.09.12-09:55:00MDT"
authorNames = "Erin L Landguth et al."

# ---------------
# Global symbols
#----------------
# when set True, routes session log traffic to BOTH the
# screen and to the log file. When False, log traffic just
# sent to log file alone.
msgVerbose = True
# File absolute paths for importing functions
SRC_PATH =  "../src/"

# ------------------------------------------
# Import Modules with Except/Try statements
# ------------------------------------------
# Python specific functions
import datetime,time,pdb,os,sys,shutil

# Numpy functions
try:
	import numpy as np                    
except ImportError as eMsg:
	print(("ImportError (%s) Numpy required."%(eMsg)))
	sys.exit(-1)

#Import the package specific folders
CDPOP_folder = os.path.dirname(os.path.abspath(SRC_PATH+"CDPOP"))

if CDPOP_folder not in sys.path:
     sys.path.insert(0, CDPOP_folder)

# CDPOP functions
try:
	from CDPOP_Modules import * 
except ImportError:
	raise ImportError("CDPOP_Modules required.")
try:
	from CDPOP_PostProcess import *
except ImportError:
	raise ImportError("CDPOP_PostProcess required.")
try:
	from CDPOP_PreProcess import *
except ImportError:
	raise ImportError("CDPOP_PreProcess required.")
try:
	from CDPOP_Mate import *
except ImportError:
	raise ImportError("CDPOP_Mate required.")
try:
	from CDPOP_Offspring import *
except ImportError:
	raise ImportError("CDPOP_Offspring required.")
try:
	from CDPOP_Disperse import *
except ImportError:
	raise ImportError("CDPOP_Disperse required.")	

#------------------------------------------------------------
# Begin main file execution
#------------------------------------------------------------ 
if __name__ == '__main__':
		
	# ------------------------------------------------------	
	# Start timer, get script arguments, create log writeout
	# ------------------------------------------------------
	# Timing events: start
	start_time = datetime.datetime.now()
	foldertime = int(time.time())
	
	if len(sys.argv) >= 4:
		datadir = sys.argv[1]+'/'
		fileans = datadir+sys.argv[2]
		outdir = datadir+sys.argv[3]+str(foldertime)+'/'
	
	# If user did not specify .rip file
	else:
		print("User must specify data directory, input file name, and output file directory (e.g., at command line type CDPOP.py ../CDPOP_data/ inputvariables16pnts.csv exampleout).")
		sys.exit(-1)	
	
	# If .ip file does not exist
	if not os.path.exists(fileans):
		print(("Cannot find or open runtime inputs file(%s)"%(fileans)))
		sys.exit(-1)
	
	# Create output file directory - will automatically put in the data directory
	os.mkdir(outdir)
	
	# This properly names log file
	logSessionPath = outdir+"cdpop.log"
	logfHndl =open(logSessionPath,'w')
	
	msgVerbose = True
	logMsg(logfHndl,"\n%s Release %s Version %s\n"%(appName,appRele,appVers))
	logMsg(logfHndl,"Author(s): %s"%(authorNames)+'\n')
	logMsg(logfHndl,"Session runtime inputs from: %s"%(fileans)+'\n\n')    
	msgVerbose = False
	
	# ------------------------------------	
	# Call DoUserInput()
	# ------------------------------------
	# Timing events: start
	start_time1 = datetime.datetime.now()
	
	# Call function and store inputvariables
	batchVars,batchVarsIndex,nSimulations = loadFile(fileans,1,',',True)
	
	# Print to log
	stringout = 'DoUserInput(): '+str(datetime.datetime.now() -start_time1) + ''
	logMsg(logfHndl,stringout)
	print('DoUserInput(): ',str(datetime.datetime.now() -start_time1),'')

	# -------------------------------------	
	# Begin Batch Looping
	# -------------------------------------
	# This loop is defined by the number of rows in inputvariables.csv
	for ibatch in range(nSimulations):
	
		# Timing events: start
		start_timeB = datetime.datetime.now()
		
		# Store all information and the type of each, also do some error checks 
		xyfilename = batchVars['xyfilename'][ibatch]
		allefreqfilename = batchVars['allefreqfilename'][ibatch]
		agefilename = batchVars['agefilename'][ibatch]
		matecdmatfile = batchVars['matecdmat'][ibatch]
		dispcdmatfile = batchVars['dispcdmat'][ibatch]
		mcruns = int(batchVars['mcruns'][ibatch])
		looptime = int(batchVars['looptime'][ibatch])
		nthfile_out = batchVars['output_years'][ibatch]
		cdclimgentimelist = batchVars['cdclimgentime'][ibatch]
		unicor_out = batchVars['output_unicor'][ibatch]		
		matemoveno = batchVars['matemoveno'][ibatch]
		matemoveparA = batchVars['matemoveparA'][ibatch]
		matemoveparB = batchVars['matemoveparB'][ibatch]
		matemoveparC = batchVars['matemoveparC'][ibatch]
		matemovethresh = batchVars['matemovethresh'][ibatch]
		freplace = batchVars['Freplace'][ibatch]
		mreplace = batchVars['Mreplace'][ibatch]
		mpaternity = batchVars['multiple_paternity'][ibatch]
		selfans = batchVars['selfans'][ibatch]
		matefreq = float(batchVars['mateFrequency'][ibatch])
		sexans = batchVars['sexans'][ibatch]
		Fdispmoveno = batchVars['Fdispmoveno'][ibatch]
		FdispmoveparA = batchVars['FdispmoveparA'][ibatch]
		FdispmoveparB = batchVars['FdispmoveparB'][ibatch]
		FdispmoveparC = batchVars['FdispmoveparC'][ibatch]
		Fdispmovethresh = batchVars['Fdispmovethresh'][ibatch]
		Mdispmoveno = batchVars['Mdispmoveno'][ibatch]
		MdispmoveparA = batchVars['MdispmoveparA'][ibatch]
		MdispmoveparB = batchVars['MdispmoveparB'][ibatch]
		MdispmoveparC = batchVars['MdispmoveparC'][ibatch]
		Mdispmovethresh = batchVars['Mdispmovethresh'][ibatch]
		philopatry = batchVars['philopatry'][ibatch]
		offnovals = batchVars['offno'][ibatch]
		Femalepercent = int(batchVars['Femalepercent'][ibatch])
		equalsexratio = batchVars['EqualsexratioBirth'][ibatch]
		twinning_pass = batchVars['TwinningPercent'][ibatch]
		popmodel = batchVars['popModel'][ibatch]
		K_envvals = batchVars['K_env'][ibatch]
		#subpopmort = batchVars['subpopmortperc'][ibatch]
		gridformat = batchVars['gridformat'][ibatch]
		muterate = float(batchVars['muterate'][ibatch])
		mutationans = batchVars['mutationtype'][ibatch]
		loci = int(batchVars['loci'][ibatch])
		intgenesans = batchVars['intgenesans'][ibatch]
		alleles = batchVars['alleles'][ibatch]
		mtdna = batchVars['mtdna'][ibatch]
		geneswap = int(batchVars['startGenes'][ibatch]) 
		cdevolveans = batchVars['cdevolveans'][ibatch]
		startSelection = int(batchVars['startSelection'][ibatch])
		betaFile_selection = batchVars['betaFile_selection'][ibatch]
		epigeneans = batchVars['epigeneans'][ibatch]
		startEpigene = int(batchVars['startEpigene'][ibatch])
		betaFile_epigene = batchVars['betaFile_epigene'][ibatch]
		cdinfect = batchVars['cdinfect'][ibatch]
		transmissionprob = float(batchVars['transmissionprob'][ibatch])
		matedist_out = batchVars['output_matedistance'][ibatch]
						
		# Distill and some error checking
		# -------------------------------
		# Grab the nthfile list range specific to user input, list or sequence
		if not isinstance(nthfile_out, (list,tuple)):
			nthfile_out = int(nthfile_out)
			if nthfile_out != 0:
				nthfile = list(range(0,looptime+nthfile_out,nthfile_out))
				del(nthfile[-1]) # Delete the last value 0, looptime - 1
			else:
				nthfile = [0]
		# If specified years with |
		else:
			nthfile = []
			# Split up list, removing space values, and appending to nthfile
			for inum in range(len(nthfile_out)):
				# Error check here if | at the end
				if len(nthfile_out[inum]) != 0:
					nthfile.append(int(nthfile_out[inum]))
		# Error check on nthfile, must be 1 less than looptime for indexing
		if max(nthfile) >= looptime:
			print('nthfile selection maximum value must be less than to looptime.')
			sys.exit(-1)
		
		# Store cdmat file information - header file (loadFile()) passes tuple or string if only 1
		if not isinstance(cdclimgentimelist, (list,tuple)):
			cdclimgentime = [cdclimgentimelist]
		else: 
			cdclimgentime = cdclimgentimelist
		if cdclimgentime[0] != '0':
			print('First cdclimate time must be 0.')
			sys.exit(-1)
	
		# Get mortality here: if tuple not returned and just one number applied across all adult ages 
		if not isinstance(agefilename,(list,tuple)):
			agefilename = [agefilename]
		else:
			agefilename = agefilename			
		if intgenesans == 'file' and allefreqfilename == 'N':
			print('Allele frequency file option specified, must give name of file.')
			sys.exit(-1)
		elif intgenesans == 'file_var' and allefreqfilename == 'N':
			print('Allele frequency file option specified, must give name of file.')
			sys.exit(-1)
		elif intgenesans == 'file' and allefreqfilename != 'N':
			if not isinstance(allefreqfilename,(list,tuple)):
				allefreqfilename = [allefreqfilename]
			else:
				allefreqfilename = allefreqfilename
		elif intgenesans == 'file_var' and allefreqfilename != 'N':
			if not isinstance(allefreqfilename,(list,tuple)):
				allefreqfilename = [allefreqfilename]
			else:
				allefreqfilename = allefreqfilename
		
		# If multiple XY files were specified for introducing individuals, then put in list as above
		if not isinstance(xyfilename,(list,tuple)):
			xyfilename = [xyfilename]
		else:	
			xyfilename = xyfilename
					
		# Create allele array
		if len(alleles.split(';')) == 1:
			alleles = int(batchVars['alleles'][ibatch])*np.ones(loci,int)
		else:
			alleles = np.asarray(alleles.split(';'),dtype = int)			
		
		# ---------------------------------
		# Some more Error checking
		# ---------------------------------
		# Have to have at least 2 alleles
		if len(np.where(alleles == 0)[0]) != 0:
			print('Must have at least 2 alleles per locus.')
			sys.exit(-1)
		if len(np.where(alleles == 1)[0]) != 0:
			print('Must have at least 2 alleles per locus.')
			sys.exit(-1)		
		
		# If cdevolve is turned on must have 2 alleles
		if cdevolveans != 'N' and alleles[0] != 2:
			print('Warning: More than 2 alleles per locus specified. CDEVOLVE only considers first 2 alleles in selection and epigenetic models, unless multiple loci models were specified.')
		if epigeneans != 'N' and alleles[0] != 2:
			print('Input Error: More than 2 alleles per locus specified. Epigenetics only considers 2 alleles in selection and epigenetic models.')
			sys.exit(-1)
		# Must have more than 1 loci
		if loci <= 1:
			print('Currently, CDPOP needs more than 1 locus to run.')
			sys.exit(-1)
			
		# Error check on forward mutation in A and backward mutation in B
		#	Can only happen if cdevolve == 2.
		if mutationans == 'forwardAbackwardBrandomN' and cdevolveans != '2':
			print('This special case of mutation is for AAbb ancestors and 2-locus selection.')
			sys.exit(-1)		
		
		# Check on parameters: equal sex ratio
		if (equalsexratio == 'N' or equalsexratio == 'AtBirth' or equalsexratio == 'WrightFisher') == False:
			print('Equal sex ratio parameter must be N, AtBirth, WrightFisher.')
			sys.exit(-1)
			
		# For female philopatry
		if philopatry == 'F' or philopatry == 'female' or philopatry == 'f' or philopatry == '0' or philopatry == 'Female': 
			
			if equalsexratio != 'AtBirth':
				print('Warning: Female philopatry is turned on and equal sex ratio at birth is recommended.')
			philopatry == 'F'
		# For male philopatry
		elif philopatry == 'M' or philopatry == 'male' or philopatry == 'm' or philopatry == '1' or philopatry == 'Male': 
			if equalsexratio != 'AtBirth':
				print('Warning: Male philopatry is turned on and equal sex ratio at birth is recommended.')
			philopatry == 'M'
		# Error if something else but N
		elif philopatry != 'N':
			print('Philopatry answer either has to be Male biased (M), Female biased (F), or unbiased (N).')
			sys.exit(-1)
					
		# grid format
		if (gridformat == 'cdpop' or gridformat == 'general' or gridformat == 'genalex' or gridformat == 'genepop' or gridformat == 'structure') == False:
			print('Grid format parameter not an option.')
			sys.exit(-1)
		
		# If genepop, some conditions
		if gridformat == 'genepop' and (len(np.where(alleles >= 99)[0]) > 0 or loci > 99):
			print('GENEPOP format requires less than 99 alleles and 99 loci.')
			sys.exit(-1)
			
		# For multiple paternity
		if mpaternity == 'Y' and (freplace != 'Y' or mreplace != 'Y'):
			print('Multiple paternity option is selected, then female and male with replacement must be both Y')
			sys.exit(-1)
		
		# Check burn in times
		if cdevolveans != 'N' and startSelection < geneswap:
			print('Start selection time must be less than genetic exchange start time (startGenes < startSelection).')
			sys.exit(-1)
		if epigeneans != 'N' and startEpigene < geneswap:
			print('Start epigenetics time must be less than genetic exchange start time (startGenes < startEpigene).')
			sys.exit(-1)
			
		# Check multiple selection model and number of loci
		if cdevolveans.split('_')[0] == 'M':
			if int(cdevolveans.split('_')[2].split('L')[1]) > loci:
				print('More loci under selection than specified number of total loci.')
				sys.exit(-1)
			if len(cdevolveans.split('_')) != 5:
				print('Multilocus selection specified, must have 6 arguments, see usermanual examples.')
				sys.exit(-1)
			if alleles[0] != int(cdevolveans.split('_')[3].split('A')[1]):
				print('Multilocus selection specified, must specify number of alleles for both neutral and selection markers.')
				sys.exit(-1)
			if intgenesans == 'random_var' or intgenesans == 'file_var':
				print('Variable allele assignment per locus can not be used with multi-locus selection.')
				sys.exit(-1)
		if epigeneans != 'N':
			if int(epigeneans.split('_')[1].split('L')[1]) > loci:
				print('More loci in epigenetic model than specified number of total loci.')
				sys.exit(-1)
			if len(epigeneans.split('_')) != 5:
				print('Epigenetic module specified, must have 4 arguments, see usermanual examples.')
				sys.exit(-1)
		if epigeneans != 'N' and cdevolveans.split('_')[0] == 'M':
			if int(cdevolveans.split('_')[2].split('L')[1]) + int(epigeneans.split('_')[1].split('L')[1]) > loci:
				print('More loci in epigenetic and selection model than specified number of total loci.')
				sys.exit(-1)		
		# Multiple files listed must equal cdclimgentime
		if len(xyfilename) > 1:
			if len(cdclimgentime) != len(xyfilename):
				print('Multiple xyfiles given, then must match the number of CDClimate generations.')
				sys.exit(-1)
			if intgenesans != 'file_introduce' and intgenesans != 'file_introduce_var':
				print('Multiple xyfiles given, then must use option file_introduce for genes.')
				sys.exit(-1)
			if len(allefreqfilename) != len(xyfilename):
				print('Multiple xyfiles given, then must match the number of allele frequency files given.')
				sys.exit(-1)
			if equalsexratio != 'N':
				print('Multiple xyfiles given, set equal sex ratio to N')
				sys.exit(-1)
			
		# Error checking for intgenesans special case for file_introduce
		if intgenesans == 'file_introduce' or intgenesans == 'file_introduce_var':
			if geneswap != 0:
				print('Gene swap starts immediately.')
		
		if popmodel == 'logistic':
			print('Logistic disabled temporarily, email erin.landguth@mso.umt.edu.')
			sys.exit(-1)
			
		# For Hindex answer
		if cdevolveans.split('_')[0] == 'Hindex':
			# Split for Gaussian
			if cdevolveans.split('_')[1] == 'Linear':
				if len(cdevolveans.split('_')[2].split(';')) != 6:
					print('CDEVOLVE answer is Hindex and 6 parameters for the Linear function must be specified, see user manual and example files.')
					sys.exit(-1)
			else:
				print('CDEVOLVE answer Hindex and only Linear option is allowed for now.')
				sys.exit(-1)
		
		# ---------------------------------------------	
		# Begin Monte-Carlo Looping
		# ---------------------------------------------
		
		# xrange(mcruns) is typically 10 - 50...and it takes a long time.
		for ithmcrun in range(mcruns):	
		
			# Timing events: start
			start_timeMC = datetime.datetime.now()
		
			# -----------------------------------------
			# Create storage variables
			# ------------------------------------------	
			# These variables will be stored in output.csv at the end of the simulation
			Population = []
			Population_age = []
			Females_age = []
			Males_age = []
			Migrants = []
			Open = []
			Track_MDeaths = []
			Track_FDeaths = []
			Births = []
			Track_MOffDeaths = []
			Track_FOffDeaths = []
			DisperseDeaths = []
			CouldNotDisperse = []
			Opt3SelectionDeaths = []
			ToTFemales = []
			ToTMales = []
			BreedFemales = []
			BreedFemales_age = []
			BreedMales = []
			Female_BreedEvents = []
			females_nomate = []
			males_nomate = []
			Alleles = []
			He = []
			Ho = []
			AllelesMutated = []
			MateDistED = []
			FDispDistED = []
			MDispDistED = []
			MateDistCD = []
			FDispDistCD = []
			MDispDistCD = []
			MateDistEDstd = []
			FDispDistEDstd = []
			MDispDistEDstd = []
			MateDistCDstd = []
			FDispDistCDstd = []
			MDispDistCDstd = []
			Infected = []
			p1 = []
			p2 = []
			q1 = []
			q2 = []
			subpopmigration = []
			subpopemigration = []
			FAvgMate = []
			MAvgMate = []
			FSDMate = []
			MSDMate = []
			MateDistances = []
			Twins = []
			Track_EpigeneMod1 = []
			Track_EpigeneMod2 = []
			Track_EpigeneDeaths = []
			Track_EpigeneReset1 = []
			Track_EpigeneReset2 = []
			maxfit = []
			minfit = []
					
			# ------------------------------------	
			# Call DoPreProcess()
			# ------------------------------------

			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			# Call function
			tupPreProcess = DoPreProcess(outdir,ibatch,ithmcrun,\
			xyfilename,agefilename,equalsexratio,loci,intgenesans,allefreqfilename,alleles,0,logfHndl,cdevolveans,cdinfect,Infected,\
			subpopmigration,subpopemigration,datadir,geneswap,epigeneans,unicor_out)
						
			ithmcrundir = tupPreProcess[0]	
			FID = tupPreProcess[1]
			id = tupPreProcess[2]
			sex = tupPreProcess[3]
			age = tupPreProcess[4]
			xgrid = tupPreProcess[5]
			xgridcopy = copy.deepcopy(xgrid)
			ygrid = tupPreProcess[6]
			ygridcopy = copy.deepcopy(ygrid)
			genes = tupPreProcess[7]
			nogrids = tupPreProcess[8]
			subpop = tupPreProcess[9]
			fitvals_pass = tupPreProcess[10]
			infection = tupPreProcess[11]
			Infected = tupPreProcess[12]
			subpopmigration = tupPreProcess[13]
			subpopemigration = tupPreProcess[14]
			Magemortvals =  tupPreProcess[15]
			Fagemortvals =  tupPreProcess[16]
			egg_lmbdavals = tupPreProcess[17]
			egg_sigmavals = tupPreProcess[18]
			allelst = tupPreProcess[19]
			Mnewmortperc = tupPreProcess[20]
			Fnewmortperc = tupPreProcess[21]
			Mmaturevals = tupPreProcess[22]
			Fmaturevals = tupPreProcess[23]
			intgenesans = tupPreProcess[24] 
			xvars_betas_pass = tupPreProcess[25]
			epimod_pass = tupPreProcess[26]
			epireset_pass = tupPreProcess[27]
			hindex = tupPreProcess[28]
			gridmort_pass = tupPreProcess[29]
						
			# Print to log
			stringout = 'DoPreProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print('DoPreProcess(): ',str(datetime.datetime.now() -start_time1),'')
			
			# -------------------------------------------
			# Start Generation Looping 
			# -------------------------------------------
			
			# Begin generation loop
			for gen in range(looptime):
			
				# Timing events: start
				start_timeGen = datetime.datetime.now()
				
				# ---------------------------------
				# Call CDClimate()
				# ---------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				# Check gen time equal to cdclimgentime
				if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:
					tupClimate = DoCDClimate(datadir,np.where(np.asarray(cdclimgentime) == str(gen))[0][0],cdclimgentime,matecdmatfile,dispcdmatfile,matemoveno,Fdispmoveno,Mdispmoveno,matemovethresh,Fdispmovethresh,Mdispmovethresh,matemoveparA,matemoveparB,matemoveparC,FdispmoveparA,FdispmoveparB,FdispmoveparC,MdispmoveparA,MdispmoveparB,MdispmoveparC,subpop,Magemortvals,Fagemortvals,offnovals,egg_lmbdavals,egg_sigmavals,K_envvals,Mnewmortperc,Fnewmortperc,fitvals_pass,twinning_pass,Mmaturevals,Fmaturevals,betaFile_selection,xvars_betas_pass,epimod_pass,epireset_pass,betaFile_epigene,cdevolveans,epigeneans,gridmort_pass)
					
					cdmatrix_mate = tupClimate[0]
					cdmatrix_F = tupClimate[1]
					cdmatrix_M = tupClimate[2]				
					thresh_mate = tupClimate[3]
					thresh_F = tupClimate[4]
					thresh_M = tupClimate[5]				
					Fdisp_ScaleMin = tupClimate[6]
					Fdisp_ScaleMax = tupClimate[7]
					Mdisp_ScaleMin = tupClimate[8]
					Mdisp_ScaleMax = tupClimate[9]
					mate_ScaleMin = tupClimate[10]
					mate_ScaleMax = tupClimate[11]
					Magemort = tupClimate[12]
					Fagemort = tupClimate[13]
					offno = tupClimate[14]
					eggs_lambda = tupClimate[15]
					eggs_sigma = tupClimate[16]
					K_env = tupClimate[17]
					Mnewmort = tupClimate[18]
					Fnewmort = tupClimate[19]
					fitvals = tupClimate[20]
					mateno = tupClimate[21]
					Fdispno = tupClimate[22]
					Mdispno = tupClimate[23]
					twinning = tupClimate[24]
					Mmature = tupClimate[25]
					Fmature = tupClimate[26]
					betas_selection = tupClimate[27]
					xvars_betas = tupClimate[28]
					epimod = tupClimate[29]
					epireset = tupClimate[30]
					betas_epigene = tupClimate[31]
					gridmort = tupClimate[32]
										
					# Error check for if nofiles == nogrids system exit
					if nogrids != len(cdmatrix_mate):
						print('The cost distance matrix dimensions are not the same as the number of individuals.')
						sys.exit(-1)
					
					# Print to log
					stringout = 'DoCDCliamte(): '+str(datetime.datetime.now() -start_time1) + ''
					logMsg(logfHndl,stringout)
					print('DoCDClimate(): ',str(datetime.datetime.now() -start_time1),'')	
					
				# -------------------------------	
				# Call ReadGrid()
				# -------------------------------
				# Use information generated from PreProcess step for first 
				#	generation, else use the following updated grid information
				if gen != 0:
					
					# Timing events: start
					start_time1 = datetime.datetime.now()
					'''move to before disperse
					# ---------------------------------
					# Add individuals if specified
					# ---------------------------------					
					if len(xyfilename) > 1:
						# Check gen time equal to cdclimgentime
						if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:
							
							# Timing events: start
							start_time1 = datetime.datetime.now()
							
							tupAddInds = AddIndividuals(cdclimgentime,gen,idnew,agenew,genesnew,sexnew,subpopnew,infectionnew,allelst,xyfilename,datadir,alleles,hindexnew)
														
							idnew = tupAddInds[0]
							sexnew = tupAddInds[1]
							agenew = tupAddInds[2]
							genesnew = tupAddInds[3]
							infectionnew = tupAddInds[4]
							hindexnew = tupAddInds[5]
							
							# Print to log
							stringout = 'AddIndividuals(): '+str(datetime.datetime.now() -start_time1) + ''
							logMsg(logfHndl,stringout)
							print(('AddIndividuals()',str(datetime.datetime.now() -start_time1),''))		
					'''
					
					tupReadGrid = ReadGrid(FIDnew,idnew,agenew,xgridnew,\
					ygridnew,genesnew,equalsexratio,sexnew,subpopnew,\
					infectionnew,allelst,geneswap,gen,intgenesans,hindexnew)	

					FID = tupReadGrid[0]
					sex = tupReadGrid[1]
					id = tupReadGrid[2]
					age = tupReadGrid[3]
					xgrid = tupReadGrid[4]
					xgridcopy = tupReadGrid[5]
					ygrid = tupReadGrid[6]
					ygridcopy = tupReadGrid[7]
					genes = tupReadGrid[8]
					nogrids = tupReadGrid[9]
					subpop = tupReadGrid[10]
					infection = tupReadGrid[11]
					filledgrids = tupReadGrid[12]					
					hindex = tupReadGrid[13]
					
					# Exit system if population is 0 or 1
					if filledgrids == 0 or filledgrids == 1:
						stringout = 'Population went extinct, program ended.'
						logMsg(logfHndl,stringout)
						print(('Population went extinct after generation '+str(gen-1)+'.\n'))
						break
					# Here we system exit if there are only F or M left
					if sexans == 'Y' and len(np.where(np.asarray(sex)=='0')[0]) == 0:
						# Print to log
						stringout = 'No females left in population, program ended.'
						logMsg(logfHndl,stringout)
						print(('There are no more females left in population after generation '+str(gen-1)+'.\n'))
						break
					if sexans == 'Y' and len(np.where(np.asarray(sex)=='1')[0]) == 0:
						# Print to log
						stringout = 'No males left in population, program ended.'
						logMsg(logfHndl,stringout)
						print(('There are no more males left in population after generation '+str(gen-1)+'.\n'))
						break				
					# Print to log
					stringout = 'ReadGrid(): '+str(datetime.datetime.now() -start_time1) + ''
					logMsg(logfHndl,stringout)
					print('ReadGrid(): ',str(datetime.datetime.now() -start_time1),'')
								
				# ---------------------------------
				# Call GetMetrics()
				# ---------------------------------
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupGetMetrics = GetMetrics(Population,nogrids,loci,alleles,genes,\
				gen,Ho,Alleles,He,subpop,p1,p2,q1,q2,Population_age,Females_age,Males_age,age,sex,Magemort,geneswap,cdevolveans,xvars_betas,betas_selection,maxfit,minfit)
				
				filledgrids = tupGetMetrics[0]
				subgridtotal = tupGetMetrics[1]
				
				# Print to log
				stringout = 'GetMetrics(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('GetMetrics(): ',str(datetime.datetime.now() -start_time1),'')
				#pdb.set_trace()
				# --------------------------------
				# Call DoEpigenetics
				# --------------------------------
				# Add time step to each tracker variable here
				if epigeneans != 'N':
					# Timing events: start
					start_time1 = datetime.datetime.now()
					
					tupEpigene = DoEpigenetics(epimod,betas_epigene,sex,id,age,genes,infection,Track_EpigeneMod1,Track_EpigeneMod2,Track_EpigeneDeaths,gen,cdevolveans,epigeneans,startEpigene,geneswap,alleles,loci)

					id = tupEpigene[0]
					sex = tupEpigene[1]
					age = tupEpigene[2]
					genes = tupEpigene[3]
					infection = tupEpigene[4]
					
					# Print to log
					stringout = 'DoEpigenetics(): '+str(datetime.datetime.now() -start_time1) + ''
					logMsg(logfHndl,stringout)
					print('DoEpigenetics(): ',str(datetime.datetime.now() -start_time1),'')
									
				# ---------------------------------------
				# Call DoMate()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupMate = DoMate(nogrids,sex,age,\
				freplace,mreplace,mateno,thresh_mate,\
				cdmatrix_mate,MateDistED,MateDistCD,xgridcopy,\
				ygridcopy,ToTMales,ToTFemales,BreedMales,BreedFemales,\
				sexans,selfans,\
				MateDistEDstd, MateDistCDstd,FAvgMate,MAvgMate,\
				FSDMate,MSDMate,filledgrids,Female_BreedEvents,gen,subpop,BreedFemales_age,Magemort,Mmature,Fmature,mate_ScaleMax,mate_ScaleMin,matemoveparA,matemoveparB,matemoveparC,MateDistances,matefreq)
				Bearpairs = tupMate[0]
				females = tupMate[1]
				females_nomate.append(tupMate[2])
				males = tupMate[3]
				males_nomate.append(tupMate[4])
				mature = tupMate[5]
				
				# Temporary, get rid of this eventually...
				CDpairs = copy.deepcopy(Bearpairs)				
				
				# Print to log
				stringout = 'DoMate(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('DoMate(): ',str(datetime.datetime.now() -start_time1),'')
						
				# ---------------------------------------
				# Call DoOffspring()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupDoOff = DoOffspring(offno,eggs_lambda,Bearpairs,CDpairs,Femalepercent,\
				Births,infection,transmissionprob,equalsexratio,\
				Mnewmort,Fnewmort,Track_MOffDeaths,Track_FOffDeaths,eggs_sigma,age,sex,twinning,Twins,subpop)
				
				offspring = tupDoOff[0]	
				offspringno = tupDoOff[1]
				
				# Print to log
				stringout = 'DoOffspring(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('DoOffspring(): ',str(datetime.datetime.now() -start_time1),'')
									
				# ---------------------------------------
				# Call InheritGenes()
				# ---------------------------------------
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				offspring = InheritGenes(gen,AllelesMutated,offspringno,offspring,genes,loci,muterate,mtdna,\
                             mutationans,geneswap,epireset,Track_EpigeneReset1,Track_EpigeneReset2,\
                             startEpigene,epigeneans,cdevolveans,alleles,hindex)
				
				# Print to log
				stringout = 'InheritGenes(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('InheritGenes(): ',str(datetime.datetime.now() -start_time1),'')
								
				# ------------------------------------------
				# Call DoAdultMortality()
				# ------------------------------------------
				#pdb.set_trace()	
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupAMort = DoMortality(nogrids,sex,id,age,xgrid,ygrid,gen,genes,Track_MDeaths,Track_FDeaths,\
                           FID,Magemort,Fagemort,infection,geneswap,popmodel,K_env,fitvals,\
                           mature,cdevolveans,Opt3SelectionDeaths,startSelection,subpop,hindex)
				
				freegrid = tupAMort[0]
				id = tupAMort[1]
				sex = tupAMort[2]
				age = tupAMort[3]
				xgrid = tupAMort[4]
				ygrid = tupAMort[5]
				genes = tupAMort[6]	
				FID = tupAMort[7]
				infection = tupAMort[8] # coming out str, change to int
				# Mature skipped, not used from here on out
				hindex = tupAMort[10]
							
				# Print to log
				stringout = 'DoAdultMortality(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('DoAdultMortality(): ',str(datetime.datetime.now() -start_time1),'')
				#pdb.set_trace()					
				# ----------------------------------------------
				# Add individuals if specified, starting gen = 1
				# ----------------------------------------------					
				if len(xyfilename) > 1 and gen != 0:
					# Check gen time equal to cdclimgentime
					if len(np.where(np.asarray(cdclimgentime) == str(gen))[0]) == 1:
						
						# Timing events: start
						start_time1 = datetime.datetime.now()
						
						#tupAddInds = AddIndividuals(cdclimgentime,gen,id,age,genes,sexnew,subpop,infection,allelst,xyfilename,datadir,alleles,hindex)
						tupAddInds = AddIndividuals(cdclimgentime,gen,id.tolist(),age.tolist(),genes,sex.tolist(),FID.tolist(),infection.tolist(),allelst,xyfilename,datadir,alleles,hindex.tolist(),subpop,freegrid)
													
						id = tupAddInds[0]
						sex = tupAddInds[1]
						age = tupAddInds[2]
						genes = tupAddInds[3]
						infection = tupAddInds[4]
						hindex = tupAddInds[5]
						FID = tupAddInds[6]
						freegrid = tupAddInds[7]
						
						# Print to log
						stringout = 'AddIndividuals(): '+str(datetime.datetime.now() -start_time1) + ''
						logMsg(logfHndl,stringout)
						print(('AddIndividuals()',str(datetime.datetime.now() -start_time1),''))
				#pdb.set_trace()			
				# ------------------------------------------
				# Call DoDisperse()
				# ------------------------------------------			
				
				# Timing events: start
				start_time1 = datetime.datetime.now()
				
				tupDoDisp = DoDisperse(offspringno,freegrid,offspring,Fdispno,\
				Mdispno,cdmatrix_F,cdmatrix_M,gen,\
				Migrants,Open,loci,alleles,\
				xgridcopy,ygridcopy,FDispDistED,MDispDistED,FDispDistCD,MDispDistCD,\
				logfHndl,cdevolveans,fitvals,FDispDistEDstd,MDispDistEDstd,\
				FDispDistCDstd,MDispDistCDstd,subpop,subpopmigration,DisperseDeaths,CouldNotDisperse,\
				gridmort,philopatry,females,subpopemigration,females_nomate[gen],\
                males,males_nomate[gen],startSelection,thresh_F,thresh_M,Fdisp_ScaleMax,\
                Fdisp_ScaleMin,Mdisp_ScaleMax,Mdisp_ScaleMin,FdispmoveparA,FdispmoveparB,\
                FdispmoveparC,MdispmoveparA,MdispmoveparB,MdispmoveparC,betas_selection,xvars_betas,maxfit,minfit)
				
				OffDisperseIN = tupDoDisp[0]
				opengrids = tupDoDisp[1]
				DispDistCD = tupDoDisp[2]
									
				# Print to log
				stringout = 'DoDisperse(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('DoDisperse(): ',str(datetime.datetime.now() -start_time1),'')
				
				# ------------------------------------------
				# Call DoOutput()
				# ------------------------------------------	
				#pdb.set_trace()
				# Timing events: start
				start_time1 = datetime.datetime.now()
							
				tupDoOut = DoOutput(nogrids,FID,OffDisperseIN,\
				xgridcopy,ygridcopy,gen,id,sex,age,xgrid,\
				ygrid,genes,nthfile,ithmcrundir,loci,alleles,subpop,\
				logfHndl,gridformat,infection,Infected,cdinfect,\
				opengrids,DispDistCD,geneswap,hindex,unicor_out)
				
				FIDnew = tupDoOut[0]
				idnew = tupDoOut[1]
				sexnew = tupDoOut[2]
				agenew = tupDoOut[3]
				xgridnew = tupDoOut[4]
				ygridnew = tupDoOut[5]
				genesnew = tupDoOut[6]
				subpopnew = tupDoOut[7]
				infectionnew = tupDoOut[8]	# mix of string/ints, fix			
				hindexnew = tupDoOut[9]
				
				# Print to log
				stringout = 'DoOutput(): '+str(datetime.datetime.now() -start_time1) + ''
				logMsg(logfHndl,stringout)
				print('DoOutput(): ',str(datetime.datetime.now() -start_time1),'')
				
				# Print to log
				stringout = 'End Generation Loop'+str(gen)+': '+str(datetime.datetime.now() -start_timeGen) + '\n'
				logMsg(logfHndl,stringout)
				print('End Generation Loop',str(gen),': ',str(datetime.datetime.now() -start_timeGen),'\n')
					
			# End::generation loop
						
			# ------------------------------------------
			# Call DoPostProcess()
			# ------------------------------------------
			
			# Timing events: start
			start_time1 = datetime.datetime.now()
			
			DoPostProcess(ithmcrundir,nogrids,\
			xgridcopy,ygridcopy,\
			loci,alleles,looptime,Population,ToTFemales,ToTMales,\
			BreedFemales,BreedMales,Migrants,Births,\
			Track_MDeaths,Track_FDeaths,Alleles,He,Ho,AllelesMutated,\
			MateDistED,FDispDistED,MDispDistED,MateDistCD,FDispDistCD,MDispDistCD,nthfile,\
			logfHndl,p1,p2,q1,q2,Infected,subpop,MateDistEDstd,\
			FDispDistEDstd,MDispDistEDstd,MateDistCDstd,FDispDistCDstd,MDispDistCDstd,subpopmigration,\
			FAvgMate,MAvgMate,FSDMate,MSDMate,DisperseDeaths,Open,CouldNotDisperse,\
			Female_BreedEvents,gridformat,subpopemigration,females_nomate,subgridtotal,Track_MOffDeaths,Track_FOffDeaths,Population_age,Females_age,Males_age,\
			BreedFemales_age,Opt3SelectionDeaths,MateDistances,matedist_out,Twins,Track_EpigeneMod1,Track_EpigeneMod2,Track_EpigeneDeaths,Track_EpigeneReset1,Track_EpigeneReset2)
			
			# Print to log
			stringout = 'DoPostProcess(): '+str(datetime.datetime.now() -start_time1) + ''
			logMsg(logfHndl,stringout)
			print('DoPostProcess(): ',str(datetime.datetime.now() -start_time1),'')
				
			# Print to log
			stringout = 'End Monte Carlo Loop'+str(ithmcrun)+': '+str(datetime.datetime.now() -start_timeMC) + '\n'
			logMsg(logfHndl,stringout)
			print('End Monte Carlo Loop',str(ithmcrun),': ',str(datetime.datetime.now() -start_timeMC),'\n')
						
		# End::Monte Carlo Loop
		
		# Print to log
		stringout = 'End Batch Loop'+str(ibatch)+': '+str(datetime.datetime.now() -start_timeB) + '\n'
		logMsg(logfHndl,stringout)
		print('End Batch Loop',str(ibatch),': ',str(datetime.datetime.now() -start_timeB),'\n')
		
	#End::Batch Loop
	
# End::Main Loop	
# Print to log
stringout = 'Total CDPOP Simulation Time: '+str(datetime.datetime.now() -start_time) + ''
logMsg(logfHndl,stringout)
logfHndl.close()
print('Total CDPOP Simulation Time: ',str(datetime.datetime.now() -start_time),'')