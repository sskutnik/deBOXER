#!/usr/bin/python3
import math
import numpy as np
import parse

# Lookup function for format specifiers for covariance matrix
def getFormat(nV,nC):
	# Convert from Fortran -> C-based indexing
	nV = nV-1
	nC = nC-1

	if(nV < 6 or nV > 13):
		raise(ValueError("Invalid value for nV; 6 < nV < 15; got {:d}".format(nV)))
	elif(nC < 0 or nC > 5):
		raise(ValueError("Invalid value for nC; 0 < nC < 7; got {:d}".format(nC)))

	ift = [ (80,'{:1d}'),(40,'{:2d}'),(26,'{:3d}'),\
		(20,'{:4d}'),(16,'{:5d}'), (13,'{:6d}'),\
		(11,'{:7.4f}'),(10,'{8.5f}'),(8,'{:9.2e}'),\
		(8,'{:10.3e}'), (7,'{:11.4e}'), (6,'{:12.5e}'),(5,'{:14.7e}')]
	return (ift[nV],ift[nC])

def getLines(myFile, numLines):
	lines = ''
	for i in range(0,numLines):
		lines += myFile.readline().strip()
	return lines

def getNumLines(nVal, valFmt):
	numLines = math.floor(nVal / valFmt[0])
	if(nVal % valFmt[0] > 0): numLines += 1			
	return numLines

def getFullMatrix(nRow, nCol, xvals, iCons):
	nx = 0
	i = 0
	j = -1
	iStart = 1
	iV = -1
	nSym = 0
	
	if(nCol == 0): 
		nSym = 1
		nCol = nRow
	
	nCon = len(iCons)
	nVal = len(xVals)

	cArray = np.zeros((nRow,nCol))

	# Load data from xvals into C(I,J) as directed by ICON
	for ic in range(0,nCon):
		if(iCons[ic] == 0):
			raise(ValueError("Invalid control value at index + {:d}"\
				.format(ic)))
		nLoad = abs(iCons[ic])
		iV += 1
		#print(iV,nVal,ic,nCon)
		assert(iV < nVal)
		cLoad = xVals[iV]

		for n in range(0, nLoad):
			j += 1
			if(nx % nVal == 0): 
				iStart += 1
				iV = 0

			if(j == nCol):
				j = 0
				i += 1
				assert (i <= nRow)
				if(nSym == 1): j = i

			if(iCons[ic] > 0):
				cLoad = 0
				if(i > iStart): cLoad = cArray[i-1][j]
			cArray[i][j] = cLoad
			nx += 1

			if(nSym == 0 or i == j): continue
			
			#print(i,j,iV)
			cArray[j][i] = cLoad
			nx +=1
	
	return(cArray)

def getBoxerData(fName, iType, mat, mt):
	headFmt = '{:1d}{:12}{:.21}' + 2*'{:5d}{:4d}' + 2*'{:4d}{:3d}' + 3*'{:4d}'
	tapeFile = open(fName, "r")
	tapeFile.seek(0,0) # reset position each time we open it
	xVals = []
	iCons = []

	# Locate the requested reaction
	for record in tapeFile:
		iTypeH = matH = 0
		parsed = parse.parse(headFmt,record)
		#print("R=",record)
		if(parsed != None):
			# Retrieve header record
			iTypeH,shortID,title,matH,mtH,mat1H,mt1H,nVal,nVf, \
				nCon, nCf, nRowM, nRowH, nColH = parsed
			if(iType == -1):
				if(matH == 0): matH = 1
				print(5*'{:6d}'.format(iTypeH, matH,mtH,mat1H, mt1H))
				break
				# Do printout of BOXER format data and continue to next record
			elif(iTypeH == 9):
				# End of file?
				exit()

		else:
			# Keep looking until we find a new header
			continue
		if(iTypeH == iType):
			# For group boundaries only, ignore MAT / MAT1  &  MT / MT1
			if(iType == 0):
				break
			if(matH == mat and mtH == mt):
				# Ignore mat1 & MT1 for iType = 1 or 2
				if(iType == 1 or iType == 2): break
				elif(iType == 3 or iType == 4):
					# Full match required for iType = 3 or 4
					if((mat1H == mat1) and (mt1H == mt1)): 
						break
			continue
		else:
			iTypeH = matH = mtH = -1
			Ci = Cj = 0


	# TODO: Figure out how to re-loop back into the parser where we stopped for iType == -1?
	if(iType > 0 and not (iTypeH == iType and matH == mat and mtH == mt)):
		print(iType, iTypeH, mat, matH, mt, mtH)

		err = "Requested data not found: " + \
			"iType = {0:d}  mat = {1:d}  MT = {2:d};  mat1 = {3:}  MT1 = {4:} ".format(iType, mat, mt, mat1, mt1)
		print(err)
		raise(RuntimeError(err))
		return

	valFmt, conFmt = getFormat(nVf, nCf)

	# Read in Boxer-formatted data for Cov matrix
	if(nVal > 0):
		numLines = getNumLines(nVal, valFmt)
		lines = getLines(tapeFile, numLines)
	
		parsed = parse.parse(nVal*valFmt[1],lines)
		if(parsed == None):
			print('nVF = {0:d}, valFmt = {1:}'.format(nVf, nVal*valFmt[1]))
			raise(ValueError)
		for item in parsed:
			xVals.append(item)

		if(valFmt[0] > 7):	
			xVals = [x * 10 for x in xVals]

	if(nCon > 0):
		numLines = getNumLines(nCon, conFmt)
		lines = getLines(tapeFile, numLines)
	
		parsed = parse.parse(nCon*conFmt[1],lines)	
		if(parsed == None):
			print('nCf = {0:d}, conFmt = {1:}'.format(nCf, conFmt))
			raise(ValueError)	
		for item in parsed:
			iCons.append(item)

	tapeFile.close()
	return ((nRowH, nColH), xVals, iCons)

inputName = './fetchCov.dat'
tapeName = "./tape71"

inputFile = open(inputName,'r')
covFile = open("newCov.dat", 'w')

#headFmt = '{:2d}{:12}{:.21}' + 11*'{:4d}'
inpFmt  = 3*'{:d}' + '{:}'

for rxn in inputFile:

	iType,mat,mt,rest = parse.parse(inpFmt,rxn)

	# mat1 & mt1 are optional; see if they're present
	tail = parse.search(2*'{:d}',rest)
	if( tail != None):
		# Assume the first two are mat1 & mt1; ignore anything else
		mat1, mt1 = rest.split()[:2]
		mat1 = int(mat1)
		mt1 = int(mt1)
	else:
		mat1 = mt1 = 0	
	#iType,shortID,title,mat,mt,mat1,mt1,*null,Ci,Cj = parse.parse(headFmt,rxn)

	if(iType == 0 and mat == 0): exit()

	try:
		dims, xVals, iCons = getBoxerData(tapeName, iType, mat, mt)
	except Exception as ex:
		print("Error parsing reaction list:\n\t{0:}".format(ex))
		exit()
	
	covArray = getFullMatrix(dims[0], dims[1], xVals, iCons)


	covFile.write(rxn)
	outStr = ''
	for row in covArray:
		for val in row:
			outStr += '{:11.4e} '.format(val) 
		outStr += '\n'
	covFile.write(outStr)

covFile.close()		

#900 IF (ITYPE.NE.-1) WRITE (NOUT,901) ITYPE,MAT,MT,MAT1,MT1
#901 FORMAT (/50H ***ERROR IN BOXR***CANNOT FIND ITYPE,MAT,MT,MAT1,
#1 ,5HMT1 =,5I5)
#IF (ITYPE.EQ.-1) WRITE(NTAB,190) IZERO,IZERO
#STOP
