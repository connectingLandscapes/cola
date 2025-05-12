# -------------------------------------------------------------------------------------------------
# CDPOP_PostProcess.py
# Author: Erin L Landguth
# Created: October 2010
# Description: This is the function/module file post processing.
# --------------------------------------------------------------------------------------------------

# Numpy functions
try:
	import numpy as np 
except ImportError:
	raise ImportError("Numpy required.")
import pdb,random,os,sys,glob,itertools
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

# --------------------------------------------------------------------------	
def nest(flat,levels):
    '''Turn a flat list into a nested list, with a specified number of lists per nesting level.
    Excess elements are silently ignored.'''
    return next(_nest(flat,levels))

def _nest(flat,levels):
    if levels:
        it = _nest(flat,levels[1:])
        while 1:
            yield list(itertools.islice(it,levels[0]))
    else:
        for d in flat:
            yield d
	
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
def DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,subpopnew,xgridnew,\
ygridnew,idnew,sexnew,agenew,genesnew,logfHndl,infection,AllDispDistCD,hindexnew,geneswap,unicor_out):
	'''
	DoGridOut_cdpop()
	Output grid.csv in cdpopformat	
	'''	
	
	# Create file to write info to
	outputfile = open(ithmcrundir+'grid'+str(gen+1)+'.csv','w')
	if unicor_out == 'Y' or unicor_out == True:
		outputfile_uni = open(ithmcrundir+'XY'+str(gen+1)+'.csv','w')
		outputfile_uni.write('XCOORD,YCOORD\n')
	
	# Write out the titles 
	title = ['Subpopulation,XCOORD,YCOORD,ID,sex,age,infection,DisperseCDist,hindex,']
	outputfile.write(title[0])	
			
	# Write out the loci title info
	# Loop for loci length
	for i in range(loci-1):
		# Loop for allele length
		for j in range(alleles[i]):
			outputfile.write('L'+str(i)+'A'+str(j)+',')
	# To get a return character on the end of the title
	for i in range(alleles[loci-1]-1):
		outputfile.write('L'+str(loci-1)+'A'+str(i)+',')
	outputfile.write('L'+str(loci-1)+'A'+str(alleles[loci-1]-1)+'\n')
	
	# Loop through each grid spot and output
	for i in range(nogrids):		
		outputfile.write(subpopnew[i]+',')
		outputfile.write(str(float(xgridnew[i]))+',')
		outputfile.write(str(float(ygridnew[i]))+',')
		outputfile.write(idnew[i]+',')
		outputfile.write(sexnew[i]+',')
		outputfile.write(str(agenew[i])+',')
		outputfile.write(str(infection[i])+',')
		outputfile.write(str(AllDispDistCD[i])+',')
		if gen >= geneswap:
			outputfile.write(str(hindexnew[i])+',')
		else:			
			outputfile.write('NA,')
		# Write out gene info
		for iall in range(sum(alleles)):
			outputfile.write(str(genesnew[i][iall])+',')
		outputfile.write('\n')
		if sexnew[i] == 'NA':
			if unicor_out == 'Y' or unicor_out == True:
				outputfile_uni.write(str(float(xgridnew[i]))+',')
				outputfile_uni.write(str(float(ygridnew[i]))+'\n')
												
	# Logging message
	stringout = 'The file grid'+str(gen+1)+'.csv has been created'
	logMsg(logfHndl,stringout)		
	
	# Close file
	outputfile.close()
	if unicor_out == 'Y' or unicor_out == True:
		outputfile_uni.close()
	
	# End::DoGridOut_cdpop()
	
# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_general(loci,alleles,ithmcrundir,logfHndl):
	'''
	DoGridOut_general()
	Output grid.csv in general genotype format	
	'''	
		
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in range(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in range(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in range(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in range(nogrids):
			genes_cdpop.append([])
			for jspot in range(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(9+sum(alleles[0:jspot])):int(9+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in range(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in range(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene general format. Email Erin.')
						sys.exit(-1)					
			
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/generalgrid'+outputfilename[1],'w')
				
		# Write out the titles that match general grid format
		title = ['Subpopulation','X','Y','ID','sex','age','infection']
			
		# Write out the title
		for ititle in range(len(title)):
			outputfile.write(title[ititle])
			outputfile.write(',')
			
		# Write out the loci title 
		for i in range(loci-1):
			outputfile.write('Locus'+str(i+1)+'a')
			outputfile.write(',')
			outputfile.write('Locus'+str(i+1)+'b')
			outputfile.write(',')
		# To get a return character on the end of the title
		outputfile.write('Locus'+str(loci-1+1)+'a')
		outputfile.write(',')
		outputfile.write('Locus'+str(loci-1+1)+'b')
		outputfile.write('\n')	
					
		# Loop through each grid spot and output
		for i in range(nogrids):
			
			outputfile.write(subpop_cdpop[i]+',')
			outputfile.write(str(float(x_cdpop[i]))+',')
			outputfile.write(str(float(y_cdpop[i]))+',')
			outputfile.write(str(id_cdpop[i])+',')
			outputfile.write(str(sex_cdpop[i])+',')
			outputfile.write(str(age_cdpop[i])+',')
			outputfile.write(str(infection_cdpop[i])+',')
			
			# Loop through each locus
			for ithloci in range(loci-1):
			
				# Loop through each allele spot at that locus
				for ithallele in range(2):
				
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+',')
				
			# Return charater on end
			outputfile.write(str(GenFormgenes[i][loci-1][0])+',')
			outputfile.write(str(GenFormgenes[i][loci-1][1])+'\n')
											
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/general'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('General grid format file conversion complete.')
	# End::DoGridOut_general()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,subgridtotal):
	'''
	DoGridOut_genalex()
	Output grid.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal)	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in range(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in range(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in range(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in range(nogrids):
			genes_cdpop.append([])
			for jspot in range(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(9+sum(alleles[0:jspot])):int(9+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in range(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in range(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene genalex format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/genalexgrid'+outputfilename[1],'w')
			
		# Write out the first and second line of GENALEX format
		outputfile.write(str(loci)+',')
		outputfile.write(str(len(genes_cdpop))+',')
		outputfile.write(str(subpopno)+',')
		outputfile.write(str(len(genes_cdpop))+'\n')
		outputfile.write(filename+'\n')		
			
		# Write out the third line of GENALEX format
		outputfile.write('Individual ID,Population,')
		for i in range(loci):
			outputfile.write('locus'+str(i+1)+'a,')
			outputfile.write('locus'+str(i+1)+'b,')
		# The trailing white space and XY information
		outputfile.write(',X,Y\n')
							
		# Loop through each grid spot and output
		for i in range(nogrids):
			
			outputfile.write('indiv'+str(i)+',')
			outputfile.write(str(subpop_cdpop[i])+',')
			
			# Loop through each locus
			for ithloci in range(loci):
			
				# Loop through each allele spot at that locus
				for ithallele in range(2):
				
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+',')
				
			# Write out trailing information - white space and x,y information
			outputfile.write(',')
			outputfile.write(str(x_cdpop[i]).strip('[').strip(']')+',')
			outputfile.write(str(y_cdpop[i]).strip('[').strip(']')+'\n')	
											
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/genalex'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('GENALEX grid format file conversion complete.')
	# End::DoGridOut_genalex()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl):
	'''
	DoGridOut_structure()
	Output grid.csv in structure genotype format	
	'''	
		
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in range(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in range(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in range(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in range(nogrids):
			genes_cdpop.append([])
			for jspot in range(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(9+sum(alleles[0:jspot])):int(9+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in range(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in range(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0':
						print('Something wrong in gene genalex format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/structuregrid'+outputfilename[1].strip('.csv')+'.stru','w')	
		
		# Write out the first line of structure format
		for i in range(loci):
			outputfile.write('locus'+str(i+1)+' ')
		outputfile.write('\n\n')
							
		# Loop through each grid spot and output
		for i in range(nogrids):
		
			# Loop through each allele spot at that locus
			for ithallele in range(2):
				
				outputfile.write(str(id_cdpop[i])+' ')
				outputfile.write(str(subpop_cdpop[i])+' ')
								
				# Loop through each locus
				for ithloci in range(loci):
					
					outputfile.write(str(GenFormgenes[i][ithloci][ithallele])+' ')
					
				# Return character
				outputfile.write('\n\n')
															
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/structure'+outputfilename[1]+'.stru has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	
	print('STRUCTURE grid format file conversion complete.')
	# End::DoGridOut_structure()

# ---------------------------------------------------------------------------------------------------	 
def DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,subgridtotal,subpop):
	'''
	DoGridOut_genalex()
	Output grid.csv in genalex genotype format	
	'''	
	subpopno = len(subgridtotal)	
	
	# Create a genes vector, appending loci information with alleles to it
	genes_genform = []
	for iloci in range(loci):
		locitemp = np.arange(1,alleles[0]+1,1)
		genes_genform.append(list(locitemp))
	
	# List files in directory
	datfileList = glob.glob(ithmcrundir+'/'+'grid*.csv')
			
	# Get length of files
	nodatfiles = len(datfileList)
	
	# Loop through each grid
	for igrid in range(nodatfiles):

		# Grab filename in List
		filename = datfileList[igrid]
				
		# Open file for reading
		inputfile = open(filename,'r')
					
		# Read lines from the file
		lines = inputfile.readlines()
					
		#Close the file
		inputfile.close()
					
		# Create an empty matrix to append to
		x = []
					
		# Split up each line in file by tab and append to empty matrix, x
		for l in lines:
			thisline = l.strip('\n').split(',')
			x.append(thisline)
			
		# And grab the number of individuals in file here
		nogrids = len(x)-1
		
		# And grab the rest of the information from the file
		sex_cdpop = []
		id_cdpop = []
		x_cdpop = []
		y_cdpop = []
		age_cdpop = []
		infection_cdpop = []
		subpop_cdpop = []
		for ispot in range(nogrids):
			subpop_cdpop.append(x[ispot+1][0])			
			x_cdpop.append(float(x[ispot+1][1]))
			y_cdpop.append(float(x[ispot+1][2]))
			id_cdpop.append(x[ispot+1][3])
			sex_cdpop.append(x[ispot+1][4])
			age_cdpop.append(x[ispot+1][5])
			infection_cdpop.append(x[ispot+1][6])
					
		# Store genetic information: genes[individual][locus][allele]
		genes_cdpop = []
		for ispot in range(nogrids):
			genes_cdpop.append([])
			for jspot in range(loci):
				genes_cdpop[ispot].append(x[ispot+1][int(9+sum(alleles[0:jspot])):int(9+sum(alleles[0:jspot+1]))])
		
		# Delete x variable
		del(x)
		
		# Store general format gene output
		GenFormgenes = []
		
		# Loop through each individual
		for ithind in range(nogrids):
			
			# Add gene individual spot 
			GenFormgenes.append([])
			
			# Loop through each locus
			for ithloci in range(loci):
			
				# Add gene individual spot 
				GenFormgenes[ithind].append([])
				
				# Loop through each allele spot at that locu
				for ithallele in range(alleles[ithloci]):
				
					# Check if allele spot is 1
					if genes_cdpop[ithind][ithloci][ithallele] == '1' or genes_cdpop[ithind][ithloci][ithallele] == '1\n':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if allele spot is 2
					elif genes_cdpop[ithind][ithloci][ithallele] == '2' or genes_cdpop[ithind][ithloci][ithallele] == '2\n':
					
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
						GenFormgenes[ithind][ithloci].append(genes_genform[ithloci][ithallele])
					
					# Check if NA in that spot
					elif genes_cdpop[ithind][ithloci][ithallele] == 'NA' or genes_cdpop[ithind][ithloci][ithallele] == 'NA\n':
						
						# Then store that unique allele number
						GenFormgenes[ithind][ithloci].append('NA')
						GenFormgenes[ithind][ithloci].append('NA')
						
					# Error check
					elif genes_cdpop[ithind][ithloci][ithallele] != '0' and genes_cdpop[ithind][ithloci][ithallele] != '0\n':
						
						print('Something wrong in gene genepop format. Email Erin.')
						sys.exit(-1)					
				
		# Create file to write matrix to
		outputfilename = filename.split('grid')
		outputfile = open(outputfilename[0]+'/genepopgrid'+outputfilename[1].strip('.csv')+'.gen','w')
			
		# Write out the first and second line of GENEPOP format
		outputfile.write(outputfilename[0]+'grid'+outputfilename[1]+'\n')
		for i in range(loci):
			outputfile.write('LOCUS-'+str(i+1)+'\n')
					
		# Write out the genes of each individual by population
		for ipop in range(subpopno):
			outputfile.write('POP\n')
			
			# Locate the index for subpop for case that is unordered.
			popindex = np.where(np.asarray(subpop) == str(ipop+1))[0]

			# Loop through each grid spot and output
			for i in range(len(popindex)):
				# Individual ID
				outputfile.write(id_cdpop[popindex[i]]+', ')								
				# Loop through each locus
				for ithloci in range(loci-1):
					templociname = ''
					# Loop through each allele spot at that locus
					for ithallele in range(2):
						if GenFormgenes[popindex[i]][ithloci][ithallele] != 'NA':
							# Add 100
							templociname = templociname + str(GenFormgenes[popindex[i]][ithloci][ithallele]+100)
						else:
							templociname = templociname + str(GenFormgenes[popindex[i]][ithloci][ithallele])
					outputfile.write(templociname+' ')
				# Loop through last loci alleles for return character
				templociname = ''
				for ithallele in range(2):	
					if GenFormgenes[popindex[i]][loci-1][ithallele] != 'NA':
						# Add 100
						templociname = templociname + str(GenFormgenes[popindex[i]][loci-1][ithallele]+100)
					else:
						templociname = templociname + str(GenFormgenes[popindex[i]][loci-1][ithallele])
				outputfile.write(templociname+'\n')		
															
		# Logging message
		stringout = 'The file grid'+outputfilename[0]+'/genepop'+outputfilename[1]+'.csv has been created'
		logMsg(logfHndl,stringout)		
		
		# Close file
		outputfile.close()
	print('GENEPOP grid format file conversion complete.')
	# End::DoGridOut_genepop()	
	
# ---------------------------------------------------------------------------------------------------	 
def DoOutput(nogrids,FID,OffDisperseIN,xgridcopy,ygridcopy,gen,\
id,sex,age,xgrid,ygrid,genes,nthfile,ithmcrundir,loci,alleles,subpop,\
logfHndl,gridformat,infection,Infected,cdinfect,opengrids,OffDispDistCD,geneswap,hindex,unicor_out):
	'''
	DoOutput()
	Generate .txt file of old+new+Immigration generations
	Input: ithmcrundir
	Output: ithmcrundir will have .csv files of x,y coord location values of
	cost distance dispersal of the old+new generation with gene info	
	'''	
		
	# ----------------------------------------------------------
	# Update grid: grid IDS either FID, OffDisperseIN, opengrids
	# -----------------------------------------------------------
	
	# checks
	if (len(FID) + len(OffDisperseIN) + len(opengrids)) != nogrids:
		print('Error with total grids.')
		sys.exit(-1)
	
	# Store new grid values
	FIDnew = []
	subpopnew = []
	xgridnew = []
	ygridnew = []
	idnew = []
	sexnew = []	
	agenew = []	
	genesnew = []
	infectionnew = []
	AllDispDistCD = []	
	hindexnew = []
	FID = np.asarray(FID)
	
	# Keep genes in within burn in period
	if gen < geneswap:
		genesnew = genes	
	
	# Get OffDisperseIN grid locations
	OffDisperseIN_gridlocs = np.asarray([OffDisperseIN[i][1] for i in range(len(OffDisperseIN))])
	
	# Loop through each grid spot to find the right output index...
	for jgrid in range(nogrids):
			
		# Find where jgrid is - FID (existing adults), OffDisperseIN (new offspring), or opengrids
		if len(np.where(FID == jgrid)[0]) == 1:
			# Where is this spot in FID, use thisindex for some vars, jgrid for others
			thisindex = np.where(FID == jgrid)[0][0]
			FIDnew.append(jgrid)
			idnew.append(id[thisindex])
			sexnew.append(str(sex[thisindex]))
			agenew.append(int(age[thisindex])+1) # age this individual
			xgridnew.append(xgridcopy[jgrid])
			ygridnew.append(ygridcopy[jgrid])
			# Write out genes info
			if gen >= geneswap:
				genesnew.append(genes[thisindex])
				hindexnew.append(hindex[thisindex])
			subpopnew.append(subpop[jgrid])
			infectionnew.append(infection[thisindex])
			AllDispDistCD.append('NoMove')			
						
		elif len(np.where(OffDisperseIN_gridlocs == jgrid)[0]) == 1:
			# Where is this spot in OffDisperseIN
			thisindex = np.where(OffDisperseIN_gridlocs == jgrid)[0][0]
			# Get this individuals information
			thisind = OffDisperseIN[thisindex]
			# Extract information to new vars
			FIDnew.append(int(thisind[1]))
			idnew.append(thisind[2])
			sexnew.append(str(int(thisind[0][4])))
			agenew.append(1) # Age is now 1
			xgridnew.append(xgridcopy[int(thisind[1])])
			ygridnew.append(ygridcopy[int(thisind[1])])
			subpopnew.append(subpop[jgrid])
			infectionnew.append(int(thisind[0][5]))
			# Write out gene info			
			if gen >= geneswap:
				thisgenes = thisind[0][7]
				#genesnew.append(nest(thisgenes,[loci,alleles[0]]))
				genesnew.append(thisgenes)
				hindexnew.append(float(thisind[0][8]))
			AllDispDistCD.append(OffDispDistCD[thisindex])			
						
		elif len(np.where(np.asarray(opengrids) == jgrid)[0]) == 1:
			# Where is this spot in opengrids, just use the jgrid index for full length vars
			# Extract information to new vars
			FIDnew.append(jgrid)
			idnew.append('OPEN')
			sexnew.append('NA')
			agenew.append('NA')
			xgridnew.append(xgridcopy[jgrid])
			ygridnew.append(ygridcopy[jgrid])
			subpopnew.append(subpop[jgrid])
			infectionnew.append('NA')
			# Write out gene info
			if gen >= geneswap:
				thisgenes = ['NA']*(sum(alleles))
				#genesnew.append(nest(thisgenes,[loci,alleles[0]]))
				genesnew.append(thisgenes)
				hindexnew.append('NA')
			AllDispDistCD.append('NA')
									
		else:
			pdb.set_trace()
			print('Grid location missing. DoOutput()')
			sys.exit(-1)
			
	# ------------------------------------------------------------
	# Write out text file for generations specified by nthfile
	# ------------------------------------------------------------
	
	# Check if nthfile == generation.
	for inthfile in range(len(nthfile)):				
		if gen == nthfile[inthfile]:			
			# Call DoGridOut_cdpop()
			DoGridOut_cdpop(ithmcrundir,gen,loci,alleles,nogrids,\
			subpopnew,xgridnew,ygridnew,idnew,sexnew,agenew,genesnew,\
			logfHndl,infectionnew,AllDispDistCD,hindexnew,geneswap,unicor_out)
	
	# Sum infectednew to Infected
	temp = []
	for i in range(len(infection)):
		if infection[i] == 1:
			temp.append(1)
	Infected.append(len(temp))
	del(temp)
	
	# Return variables from this argument
	tupDoOut = FIDnew,idnew,sexnew,agenew,xgridnew,ygridnew,\
	genesnew,subpopnew,infectionnew,hindexnew
	return tupDoOut
	
	# End::DoOutput()
	
# ---------------------------------------------------------------------------------------------------	 
def DoPostProcess(ithmcrundir,nogrids,\
xgridcopy,ygridcopy,\
loci,alleles,looptime,Population,ToTFemales,ToTMales,\
BreedFemales,BreedMales,Migrants,Births,\
MDeaths,FDeaths,Alleles,He,Ho,AllelesMutated,\
MateDistED,FDispDistED,MDispDistED,MateDistCD,FDispDistCD,\
MDispDistCD,nthfile,logfHndl,p1,p2,q1,q2,Infected,subpop,\
MateDistEDstd,FDispDistEDstd,MDispDistEDstd,MateDistCDstd,\
FDispDistCDstd,MDispDistCDstd,subpopmigration,FAvgMate,MAvgMate,\
FSDMate,MSDMate,DisperseDeaths,Open,CouldNotDisperse,\
Female_BreedEvents,gridformat,subpopemigration,females_nomate,\
subgridtotal,MOffDeaths,FOffDeaths,Population_age,Females_age,Males_age,BreedFemales_age,SelectionDeaths,MateDistances,matedist_out,Twins,Track_EpigeneMod1,Track_EpigeneMod2,Track_EpigeneDeaths,Track_EpigeneReset1,Track_EpigeneReset2):
	'''
	DoPostProcess()
	Create Distance Matrices - Geographic, Genetic, and Cost
	and output.csv file.
	'''	
		
	# ------------------------
	# Grid format options
	# ------------------------
	# General format
	if gridformat == 'general':		
		DoGridOut_general(loci,alleles,ithmcrundir,logfHndl)
		
	# GENALEX format
	elif gridformat == 'genalex':
		DoGridOut_genalex(loci,alleles,ithmcrundir,logfHndl,subgridtotal)
		
	# STRUCTURE format
	elif gridformat == 'structure':
		DoGridOut_structure(loci,alleles,ithmcrundir,logfHndl)
	
	# GENEPOP format
	elif gridformat == 'genepop':
		DoGridOut_genepop(loci,alleles,ithmcrundir,logfHndl,subgridtotal,subpop)
			
	# -------------------------------------------------
	# output.csv
	# -------------------------------------------------
			
	# Calculate other metrics here
	# ----------------------------
	tempPop = np.asarray(Population,dtype='float')[:,0]
	growthPop = tempPop[1:]/tempPop[0:(len(tempPop)-1)]
	
	# Create time array
	time = np.arange(0,len(Population),1)
		
	# Get unique number of subpops
	nosubpops = len(np.unique(subpop))
	
	# Store Ouput in file
	# Create file to write info to
	outputfile = open(ithmcrundir+'output.csv','w')
	
	# Write out the titles
	# Add Titles from xypoints
	if matedist_out == 'Y':
		outputtitle = ['Year','Population','Population_Age1+','GrowthRate','ToTFemales','ToTFemales_Age1+','ToTMales','ToTMales_Age1+','BreedFemales','BreedFemales_Age1+','BreedMales','BreedEvents_Females','Females_NoMate','Migrants','SelectionDeaths','Births','Male_Age0Deaths','Female_Age0Deaths','Male_AgeDeaths1+','Female_AgeDeaths1+','Alleles','He','Ho','Mutations','MateDistED','MateDistEDstd','Female_DispDistED','Female_DispDistEDstd','Male_DispDistED','Male_DispDistEDstd','MateDistCD','MateDistCDstd','Female_DispDistCD','Female_DispDistCDstd','Male_DispDistCD','Male_DispDistCDstd','p1','p2','q1','q2','Infected','SubpopImmigration','SubpopEmigration','SubpopNoMate','FemalesMeanMate','MalesMeanMate','FemalesSDMate','MalesSDMate','OpenLocations','CouldNotDisperse','MatureSelectionDeaths','Twins','EpigeneMod_A1','EpigeneMod_A2','EpigeneDeaths','EpigeneResets_A1','EpigeneResets_A2','AllMateCDistances']
	else:
		outputtitle = ['Year','Population','Population_Age1+','GrowthRate','ToTFemales','ToTFemales_Age1+','ToTMales','ToTMales_Age1+','BreedFemales','BreedFemales_Age1+','BreedMales','BreedEvents_Females','Females_NoMate','Migrants','DisperseDeaths','Births','Male_Age0Deaths','Female_Age0Deaths','Male_AgeDeaths1+','Female_AgeDeaths1+','Alleles','He','Ho','Mutations','MateDistED','MateDistEDstd','Female_DispDistED','Female_DispDistEDstd','Male_DispDistED','Male_DispDistEDstd','MateDistCD','MateDistCDstd','Female_DispDistCD','Female_DispDistCDstd','Male_DispDistCD','Male_DispDistCDstd','p1','p2','q1','q2','Infected','SubpopImmigration','SubpopEmigration','SubpopNoMate','FemalesMeanMate','MalesMeanMate','FemalesSDMate','MalesSDMate','OpenLocations','CouldNotDisperse','MatureSelectionDeaths','Twins','EpigeneMod_A1','EpigeneMod_A2','EpigeneDeaths','EpigeneResets_A1','EpigeneResets_A2']
		
	
	# Write out the title
	for i in range(len(outputtitle)-1):
		outputfile.write(outputtitle[i])
		outputfile.write(',')
	# To get return character on the end
	outputfile.write(str(outputtitle[len(outputtitle)-1]))				
	outputfile.write('\n')
	
	# Write to file
	for i in range(len(time)):
		outputfile.write(str(time[i])+',')
		for j in range(nosubpops+1):
			outputfile.write(str(Population[i][j])+'|')
		outputfile.write(',')
		for j in range(len(Population_age[i])):#Split by population
			for iage in range(len(Population_age[0][0])-1):#Split by age
				outputfile.write(str(Population_age[i][j][iage])+';')				
			outputfile.write(str(Population_age[i][j][len(Population_age[0][0])-1])+'|')
		outputfile.write(',')		
		if i != len(time)-1:
			outputfile.write(str(growthPop[i])+',')
		else:
			outputfile.write('NA,')
		for j in range(nosubpops+1):
			outputfile.write(str(ToTFemales[i][j])+'|')
		outputfile.write(',')
		for j in range(len(Females_age[i])):#Split by population
			for iage in range(len(Females_age[0][0])-1):#Split by age
				outputfile.write(str(Females_age[i][j][iage])+';')				
			outputfile.write(str(Females_age[i][j][len(Females_age[0][0])-1])+'|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(ToTMales[i][j])+'|')
		outputfile.write(',')
		for j in range(len(Males_age[i])):#Split by population
			for iage in range(len(Males_age[0][0])-1):#Split by age
				outputfile.write(str(Males_age[i][j][iage])+';')				
			outputfile.write(str(Males_age[i][j][len(Males_age[0][0])-1])+'|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(BreedFemales[i][j])+'|')
		outputfile.write(',')
		for j in range(len(BreedFemales_age[i])):
			outputfile.write(str(BreedFemales_age[i][j])+'|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(BreedMales[i][j])+'|')
		outputfile.write(',')		
		outputfile.write(str(Female_BreedEvents[i])+',')
		outputfile.write(str(len(females_nomate[i]))+',')
		outputfile.write(str(Migrants[i])+',')
		outputfile.write(str(DisperseDeaths[i])+',')
		outputfile.write(str(Births[i])+',')
		outputfile.write(str(MOffDeaths[i])+',')
		outputfile.write(str(FOffDeaths[i])+',')
		for j in range(len(MDeaths[i])):#Split by population
			for iage in range(len(Population_age[0][0])-1):#Split by age
				if len(MDeaths[i][j]) == len(Population_age[0][0]):
					outputfile.write(str(MDeaths[i][j][iage])+';')
				else:
					outputfile.write('NA;')
			if len(MDeaths[i][j]) == (len(Population_age[0][0])):
				outputfile.write(str(MDeaths[i][j][len(Population_age[0][0])-1])+'|')
			else:
				outputfile.write('NA|')
		outputfile.write(',')
		for j in range(len(FDeaths[i])):#Split by population
			for iage in range(len(Population_age[0][0])-1):#Split by age 
				if len(FDeaths[i][j]) == len(Population_age[0][0]):
					outputfile.write(str(FDeaths[i][j][iage])+';')
				else:
					outputfile.write('NA;')
			if len(FDeaths[i][j]) == (len(Population_age[0][0])):
				outputfile.write(str(FDeaths[i][j][len(Population_age[0][0])-1])+'|')
			else:
				outputfile.write('NA|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(Alleles[i][j])+'|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(He[i][j])+'|')
		outputfile.write(',')
		for j in range(nosubpops+1):
			outputfile.write(str(Ho[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(AllelesMutated[i])+',')
		outputfile.write(str(MateDistED[i])+',')
		outputfile.write(str(MateDistEDstd[i])+',')
		outputfile.write(str(FDispDistED[i])+',')
		outputfile.write(str(FDispDistEDstd[i])+',')
		outputfile.write(str(MDispDistED[i])+',')
		outputfile.write(str(MDispDistEDstd[i])+',')
		outputfile.write(str(MateDistCD[i])+',')
		outputfile.write(str(MateDistCDstd[i])+',')
		outputfile.write(str(FDispDistCD[i])+',')
		outputfile.write(str(FDispDistCDstd[i])+',')
		outputfile.write(str(MDispDistCD[i])+',')
		outputfile.write(str(MDispDistCDstd[i])+',')
		outputfile.write(str(p1[i])+',')
		outputfile.write(str(p2[i])+',')
		outputfile.write(str(q1[i])+',')
		outputfile.write(str(q2[i])+',')
		outputfile.write(str(Infected[i])+',')
		for j in range(nosubpops):
			outputfile.write(str(subpopmigration[i][j])+'|')
		outputfile.write(',')
		for j in range(nosubpops):
			outputfile.write(str(subpopemigration[i][j])+'|')
		outputfile.write(',')
		outputfile.write(str(FAvgMate[i])+',')
		outputfile.write(str(MAvgMate[i])+',')
		outputfile.write(str(FSDMate[i])+',')
		outputfile.write(str(MSDMate[i])+',')
		outputfile.write(str(Open[i])+',')
		outputfile.write(str(CouldNotDisperse[i])+',')
		outputfile.write(str(SelectionDeaths[i])+',')
		outputfile.write(str(Twins[i])+',')
		if len(Track_EpigeneMod1) != 0:
			outputfile.write(str(Track_EpigeneMod1[i])+',')
			outputfile.write(str(Track_EpigeneMod2[i])+',')
			outputfile.write(str(Track_EpigeneDeaths[i])+',')
			outputfile.write(str(Track_EpigeneReset1[i])+',')
			outputfile.write(str(Track_EpigeneReset2[i])+',')
		else:
			outputfile.write('NA,')
			outputfile.write('NA,')
			outputfile.write('NA,')
			outputfile.write('NA,')
			outputfile.write('NA,')
		if matedist_out == 'Y':
			outputfile.write(str(MateDistances[i]))
			#for j in xrange(len(MateDistances[i])):
			#	outputfile.write(str(MateDistances[i][j])+'|')
			#outputfile.write(',')
		outputfile.write('\n')
				
	# Logging message
	stringout = 'The file outputfile.csv has been created'
	logMsg(logfHndl,stringout)	
	
	# Close file
	outputfile.close()