# ---------------------------------------------------------------------------------------------------	
def DoDaGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoDaGeneticDistance()
	This function outputs the genetic distance matrix 
	following Nei's algorithm.
	'''	
	
	print('Nei Da is not functional yet for versions greater than 1.2. Email Erin.')
	sys.exit(-1)
	
	# Insert -1 at beginning of nthfile, to get grid0.csv
	nthfile.insert(0,-1)
	
	# Get the specified nthfile's grid to calculated matrix
	for i in range(len(nthfile)):
						
		# Open file for reading
		inputfile = open(ithmcrundir+'grid'+str(nthfile[i]+1)+'.csv','r')
		
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in range(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in range(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(float(tempgenes[k][j]))
		
		# Create a matrix of zeros to be filled
		gendmatrix = np.zeros((nogrids,nogrids),float)
		
		# Loop through each individual k
		for k in range(nogrids):
			# Compare individual k to every other inidividual j
			for j in range(nogrids):
				# Create a tempvariable to be written over for each comparison
				tempsqrt=[]
				# Loop through each locus
				for locus in range(loci):						
					# Loop through each allele value
					for alle in range(alleles[locus]):
						# Find the allele frequencies between k and j checking the 4 conditions
						if genes[k][alle+sum(alleles[0:locus])]==2.0:
							if genes[j][alle+sum(alleles[0:locus])]==2.0:
								tempsqrt.append(np.sqrt(1.0))
							elif genes[j][alle+sum(alleles[0:locus])]==1.0:
								tempsqrt.append(np.sqrt(1.0*0.5))
						elif genes[k][alle+sum(alleles[0:locus])]==1.0:
							if genes[j][alle+sum(alleles[0:locus])]==1.0:
								tempsqrt.append(np.sqrt(0.5*0.5))
							elif genes[j][alle+sum(alleles[0:locus])]==2.0:
								tempsqrt.append(np.sqrt(0.5*1.0))
				# Write the Da value to gendmatrix
				gendmatrix[k][j] = 1-float(sum(tempsqrt))/(loci)
		
		# Transpose cdmatrix to get vertical display (just being picky)
		gendmatrix = np.transpose(gendmatrix)
		
		# Create file to write matrix to
		outputfile = open(ithmcrundir+'Gdmatrix'+str(nthfile[i]+1)+'.csv','w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:
		
			# Grab each element in each row and write element to outputfile
			for ele in range(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')
		
		# Close file
		outputfile.close()
		
		# Logging message
		stringout = 'Gdmatrix'+str(nthfile[i])+'.csv has been created'
		logMsg(logfHndl,stringout)
	
	# Delete the -1
	del(nthfile[0])
	
	# End::DoDaGeneticDistance()

# ---------------------------------------------------------------------------------------------------	
def DoDpsGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoDpsGeneticDistance()
	This function outputs the genetic distance matrix 
	following proportion of shared alleles algorithm.
	'''	
	# List the grid files in directory
	csvfileList = glob.glob(ithmcrundir+'grid*.csv')
		
	# Get the specified nthfile's grid to calculated matrix
	for i in range(len(csvfileList)):
						
		# Open file for reading
		inputfile = open(csvfileList[i],'r')
			
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in range(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in range(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(tempgenes[k][j].strip('\n'))
		
		# Create a matrix of zeros to be filled
		gendmatrix = []
		
		# Loop through each individual k
		tempcount = 0		
		for k in range(nogrids):
			
			# Break if NA
			if genes[k][0] == 'NA':
				continue
			else:
				gendmatrix.append([])
				tempcount = tempcount+1
				# Compare individual k to every other inidividual j
				for j in range(nogrids):
					
					# Break if NA
					if genes[j][0] == 'NA':
						continue
					else:
						# Create a tempvariable to be written over for each comparison
						tempmin=[]
						# Loop through each allele value
						for alle in range(sum(alleles)):
							# Find the shared alleles between k and j checking the 4 conditions
							if genes[k][alle]=='2':
								if genes[j][alle]=='2':
									tempmin.append(2)
								elif genes[j][alle]=='1':
									tempmin.append(1)
							elif genes[k][alle]=='1':
								if genes[j][alle]=='1':
									tempmin.append(1)
								elif genes[j][alle]=='2':
									tempmin.append(1)
						# Write the Dps value to gendmatrix
						gendmatrix[tempcount-1].append(1-float(sum(tempmin))/(2*loci))
		
		# Strip directory/filename of grid and add 'Gdmatrix.csv'
		gdpathname = csvfileList[i].replace('grid','Gdmatrix')
	
		# Create file to write matrix to
		outputfile = open(gdpathname,'w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:
		
			# Grab each element in each row and write element to outputfile
			for ele in range(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')
		
		# Close file
		outputfile.close()
		
		# Logging message
		stringout = gdpathname+' has been created'
		logMsg(logfHndl,stringout)
		print(('The genetic distance matrix '+gdpathname+' has been created.'))
		
	# End::DoDpsGeneticDistance()
		
# ---------------------------------------------------------------------------------------------------	
def DoBrayCurtisGeneticDistance(loci,nogrids,alleles,\
ithmcrundir,nthfile,logfHndl):
	'''
	DoBrayCurtisGeneticDistance()
	This function outputs the genetic distance matrix 
	following the Bray-Curtis algorithm.
	'''	
	print('Bray-Curtis is not functional yet for versions greater than 1.2. Email Erin.')
	sys.exit(-1)
	A = loci*2
	B = loci*2
	
	# Insert -1 at beginning of nthfile, to get grid0.csv
	nthfile.insert(0,-1)
	
	# Get the specified nthfile's grid to calculated matrix
	for i in range(len(nthfile)):
						
		# Open file for reading
		inputfile = open(ithmcrundir+'grid'+str(nthfile[i]+1)+'.csv','r')
			
		# Read lines from the file
		lines = inputfile.readlines()
		
		#Close the file
		inputfile.close()
		
		# Create an empty matrix to append to
		x = []
		
		# Split up each line in file and append to empty matrix, x
		for l in lines:
			thisline = l.split(',')
			x.append(thisline)
				
		# Store genetic information: genes[individual], but need them as float values
		genes = []
		tempgenes = []
		for k in range(len(x)-1):
			# Get list from read in file
			tempgenes.append(x[k+1][8:int(8+sum(alleles))])
			# Create spot in genes
			genes.append([])
			for j in range(sum(alleles)):
				# Make each list spot an integer
				genes[k].append(float(tempgenes[k][j]))	
		
		# Create a matrix of zeros to be filled
		gendmatrix = np.zeros((nogrids,nogrids),float)

		# Loop through the rows
		for igrids in range(nogrids):
			
			# Loop through the columns
			for j in range(nogrids):
				
				# Create a tempvariable to be written over
				tempW=[]
				
				# Loop through each allele value
				for k in range(sum(alleles)):
					
					# Find the minimum value between the i and jth entry
					tempW.append(min(int(genes[igrids][k]),int(genes[j][k])))
					
				# Write the Curtis-Bray value to gendmatrix
				gendmatrix[igrids][j] = 1-(float((2*sum(tempW)))/(A+B))

		# Create file to write matrix to
		outputfile = open(ithmcrundir+'Gdmatrix'+str(nthfile[i]+1)+'.csv','w')
		
		# Sequence each row in the matrix
		for seqrow in gendmatrix:

			# Grab each element in each row and write element to outputfile
			for ele in range(len(seqrow)):
				outputfile.write(str(seqrow[ele]))
				# Add comma
				outputfile.write(',')
			
			# Return line
			outputfile.write('\n')

		# Close file
		outputfile.close()

		# Logging message
		stringout = 'Gdmatrix'+str(nthfile[i])+'.csv has been created'
		logMsg(logfHndl,stringout)
	
	# Delete the -1
	del(nthfile[0])
	
	# End::DoBrayCurtisGeneticDistance()
