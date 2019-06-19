#!/usr/bin/python3
import sys, traceback
from collections import OrderedDict
import warnings
import math
import numpy as np
import parse

def getLinesToSkip(fmt, nVC):
	valFmt, conFmt = getFormat(fmt)
	skipLines = math.ceil(nVC[0] / valFmt[0]) + \
		math.ceil(nVC[1] / conFmt[0])
	return skipLines

# Get list of requested reactions and return a dictionary
def getReactionList(inputName):

	inputFile = open(inputName,'r')
	inpFmt  = 3*'{:d}' + '{:}'

	rxnDict = {}
	for rxn in inputFile:
		iType = 0
		mat1 = 0
		mt1 = 0
		parsed = parse.parse(inpFmt,rxn)
	
		if(parsed != None):
			iType,mat,mt,rest = parsed
		else:
			parsed = rxn.split()
			if(parsed != None and len(parsed) >= 2):
				iType = int(parsed[0])
				mat = int(parsed[1])
			else:
				print("Can't parse reaction: {0:}".format(rxn))
				continue
	
		# mat1 & mt1 are optional; see if they're present
			tail = parse.search(2*'{:d}',rest)
			if( tail != None):
				# Assume the first two are mat1 & mt1; ignore anything else
				mat1, mt1 = rest.split()[:2]
				mat1 = int(mat1)
				mt1 = int(mt1)
			else:
				mat1 = mt1 = 0	
		
			if(iType == 0 and mat == 0): 
				print("Found terminator; quitting")
				exit()
	
		# Create a dictionary of materials requested, grouping multiple 
		# reactions into one common material
		if(iType == 3 or iType == 4): 
			mt1 = mt
			mat1 = mat
	
		thisRxn =  {'MT': mt, 'iType': iType, 'MAT1': mat1, 'MT1': mt1 }
		if(mat not in rxnDict):
			rxnDict[mat] = [thisRxn]
		else:
			rxnDict[mat].append(thisRxn)

	inputFile.close()
	return rxnDict

# Lookup function for format specifiers for covariance matrix
def getFormat(fmt):
	nV, nC = fmt
	# Convert from Fortran -> C-based indexing
	nV -= 1
	nC -= 1

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
		lines += (myFile.readline().strip() + ' ')

	return lines.strip()

def getNumLines(nVal, valFmt):
	numLines = math.floor(nVal / valFmt[0])
	if(nVal % valFmt[0] > 0): numLines += 1			
	return numLines

#def getFullMatrix(nRow, nCol, xvals, iCons):
def getFullMatrix(tapeFile, iType, mat, mt, mat1=0, mt1=0, iStart = 0, cArray=None):

	dims, fmt, nVC, nRowM = getBoxerHeader(tapeFile, iType, mat, mt, mat1, mt1)
	xVals, iCons = getBoxerData(tapeFile, fmt, nVC)	

	nx = 0
	i = 0
	j = -1
	iV = 0
	nSym = 0
		
	nVal = nVC[0]
	nCon = nVC[1]
	nRow, nCol = dims

	if(nCol == 0): 
		nSym = 1
		nCol = nRow

	if(cArray is None): cArray = np.zeros((nRow,nCol))

	# Load data from xvals into C(I,J) as directed by ICON
	for ic in range(0,nCon):
		if(iCons[ic] == 0):
			#print(iCons)
			raise(ValueError("Invalid control value at index + {:d}"\
				.format(ic)))


		nLoad = abs(iCons[ic])

		if(iCons[ic] < 0): 
			iV += 1
			assert(iV < nVal)
			cLoad = xVals[iV]

		for n in range(0, nLoad):
			j += 1
			if(nx % nVal == 0): 
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
			
			cArray[j][i] = cLoad
			nx += 1

	if(nRowM > 0):
		# Get next page of data
		iStart += 1	
		return getFullMatrix(tapeFile, iType, mat, mt, mat1, mt1, iStart, cArray)
	else:
		return cArray	


def getBoxerHeader(tapeFile, iType, mat, mt, mat1=0, mt1=0):
	headFmt = '{:1d}{:12}{:.21}' + 2*'{:5d}{:4d}' + 2*'{:4d}{:3d}' + 3*'{:4d}'

	# Locate the requested reaction
	for record in tapeFile:
		iTypeH = matH = 0
		parsed = parse.parse(headFmt,record)
		if(not parsed):
			# Keep looking until we find a new header
			print("Parse failed: ", record)
			continue

		iTypeH,shortID,title,matH,mtH,mat1H,mt1H,nVal,nVf, \
			nCon, nCf, nRowM, nRowH, nColH = parsed

		if(iType == -1):
			if(matH == 0): matH = 1
			print(5*'{:6d}'.format(iTypeH, matH,mtH,mat1H, mt1H))
			break
			# Do printout of BOXER format data and continue to next record
		elif(iTypeH == 9):
			# End of file?
			message = ("Reaction iType = {0:d} / mat = {1:d}" + \
				   " / MT = {2:d} not found! Continuing...") \
				  .format(iType, mat, mt)
			warnings.warn(message)
			tapeFile.close()
			return (-1, None, None)


		if(iTypeH == iType):
			# For group boundaries only, ignore MAT / MAT1  &  MT / MT1
			if(iType == 0):
				break
			if(matH == mat and mtH == mt):
				# Ignore mat1 & MT1 for iType = 1 or 2
				if(iType == 1 or iType == 2): 
					break
				elif(iType == 3 or iType == 4):
					# Full match required for iType = 3 or 4
					if((mat1H == mat1) and (mt1H == mt1)): 
						break
		else:
			iTypeH = matH = mtH = -1
			Ci = Cj = 0
		# Skip to the next header record
		toSkip = getLinesToSkip((nVf, nCf), (nVal, nCon))

		for i in range(0,toSkip): next(tapeFile)
		continue

	# TODO: Figure out how to re-loop back into the parser where we stopped for iType == -1?
	if(iType > 0 and not (iTypeH == iType and matH == mat and mtH == mt)):

		err = "Requested data not found: " + \
			"\n\tiType = {0:d}\n\tmat = {1:d}  MT = {2:d};\n\tmat1 = {3:d}  MT1 = {4:d}" \
			.format(iType, mat, mt, mat1, mt1)
		#print(err)
		raise(RuntimeError(err))
		return
	return ((nRowH, nColH), (nVf, nCf), (nVal, nCon), nRowM)


def getBoxerData(tapeFile, fmt, nVC):
	nVal, nCon = nVC

	xVals = []
	iCons = []
	valFmt, conFmt = getFormat(fmt) # fmt = (nVf, nCf)

	# Read in Boxer-formatted data for Cov matrix
	if(nVal > 0):
		numLines = getNumLines(nVal, valFmt)
		lines = getLines(tapeFile, numLines)
		parsed = parse.parse(nVal*valFmt[1],lines)
		#print(parsed)
		if(parsed is None):
			print('nVF = {0:d}, valFmt = {1:}'.format(fmt[0], nVal*valFmt[1]))
			raise(ValueError)
			return
		for item in parsed:
			xVals.append(item)

		if(valFmt[0] > 7):	
			xVals = [x * 10 for x in xVals]

	if(nCon > 0):
		numLines = getNumLines(nCon, conFmt)
		lines = getLines(tapeFile, numLines)
		#parsed = parse.parse(nCon*conFmt[1],lines)	
		parsed = lines.split()
		
		if(parsed == None):
			print('nCf = {0:d}, conFmt = {1:}'.format(fmt[1], conFmt))
			raise(ValueError)	
		for item in parsed:
			iCons.append(int(item))

	return (xVals, iCons)


tapeName = "./tape81"
inputName = './fetchCov.dat'

rxnDict = getReactionList(inputName)
covFile = open("newCov.dat", 'w')

firstRxn = True

maxRxns = max(len(rxnDict[key]) for key in rxnDict)
numMat = len(rxnDict.keys())
totRxns = sum(len(rxnDict[key]) for key in rxnDict)

rxnDict = OrderedDict(sorted(rxnDict.items()))
tapeFile = open(tapeName, "r")

for key in rxnDict.keys():

	if(firstRxn):
		try:
			dims, fmt, nVC, nRowM = getBoxerHeader(tapeFile, 0, key, 0)	
			xVals, iCons = getBoxerData(tapeFile, fmt, nVC)	

		except Exception as ex:
			print("Problem getting energy bounds for material:\t{0:}".format(key))
			exit()

		firstRxn = False
		xVals = np.array(xVals)
	
		covFile.write('\t{0:d}\tES12.5\n'.format(len(xVals)))
	
			
		bounds = ''.join('{' + '{:d}:11.4e'.format(i) + '} ' for i in range(0,len(xVals)))
		covFile.write(bounds.format(*xVals))
		covFile.write('\n\t{0:d}\t{1:d}\t{2:d}\n'.format(numMat,totRxns,maxRxns))

	covFile.write('\t{0:d}\n'.format(len(rxnDict[key])))

	for rxn in rxnDict[key]:	
		# Grab up energy bounds for first reaction & output to master file
		print("Fetching MAT = {0:d} & MT = {1:d} with iType = {2:d}".format(key, rxn['MT'],rxn['iType']))			
		tapeFile.seek(0,0) # reset position each time we open it
		
		covArray = None
		#try:
		covArray = getFullMatrix(tapeFile, rxn['iType'], key, rxn['MT'], rxn['MAT1'],rxn['MT1'], 0, covArray)
		#except Exception as ex:
		#	print("Error parsing reaction list:\n\t{0:}".format(ex))
			#traceback.print_exc(file=sys.stdout)
		#	exit()
		#if(dims == -1):
		#	continue
		assert(covArray is not None)
	
		covFile.write('\t{0:d}\t{1:d}\t{2:d}\n'.format(key,rxn['MT'],rxn['iType']))
		outStr = ''
		for row in covArray:
			for val in row:
				outStr += '{:11.4e} '.format(val) 
			outStr += '\n'

		covFile.write(outStr)

tapeFile.close()
covFile.close()		

#bounds = ''.join('{' + '{:d}:11.4e'.format(i) + '} ' for i in range(0,6))
	#for row in zip(*[iter(xVals)]*6):
	#	print(bounds.format(*row))
	#	covFile.write(bounds.format(*row) + "\n")	
	# Handle the last row (which may have < 6 elements)
	#lastBounds = len(xVals) % 6 
	#bounds = ''.join('{' + '{:d}:12.5e'.format(i) + '} ' for i in range(0,lastBounds))
	#covFile.write(bounds.format(*xVals[-lastBounds:]))
