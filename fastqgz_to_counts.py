# pipeline to generate read counts and phenotype scores directly from gzipped sequencing data

#FEATURE to add: sequential unmapped read alignment

import os
import sys
import subprocess
import gzip
import multiprocessing
import fnmatch
import glob
import csv
import argparse
from Bio import SeqIO


### Sequence File to Trimmed Fasta Functions ###

def parallelSeqFileToTrimmedFasta(fastqGzFileNameList, fastaFileNameList, outbase, processPool, startIndex=None, stopIndex=None, test=False):
	"""Function Description

	Parameters
	----------


	Returns
	-------


	Raises
	------


	"""

	if len(fastqGzFileNameList) != len(fastaFileNameList):
		raise ValueError('In and out file lists must be the same length')

	arglist = zip(fastqGzFileNameList, fastaFileNameList, [startIndex]*len(fastaFileNameList),[stopIndex]*len(fastaFileNameList), [test]*len(fastaFileNameList))
	#print arglist

	readsPerFile = processPool.map(seqFileToTrimmedFastaWrapper, arglist)

	return zip(fastaFileNameList,readsPerFile)


def seqFileToTrimmedFastaWrapper(arg):
	
	trimmingFuction = None

	for fileTup in acceptedFileTypes:
		if fnmatch.fnmatch(arg[0],fileTup[0]):
			trimmingFuction = fileTup[1]

	if trimmingFuction == None:
		raise ValueError('Sequencing file type not recognized!')

	return trimmingFuction(*arg)

def fastqGzToTrimmedFasta(fastqGzFileName, fastaFileName, startIndex=None, stopIndex=None, test=False):
	with gzip.open(fastqGzFileName) as filehandle:
		#fastq file subtype shouldn't matter since quality scores are discarded
		numReads = seqFileToTrimmedFasta_filehandle(filehandle,fastaFileName,'fastq',startIndex,stopIndex,test)

	return numReads

def fastqToTrimmedFasta(fastqFileName, fastaFileName, startIndex=None, stopIndex=None, test=False):
	with open(fastqFileName) as filehandle:
		numReads = seqFileToTrimmedFasta_filehandle(filehandle,fastaFileName,'fastq',startIndex,stopIndex,test)

	return numReads

def fastaToTrimmedFasta(inFastaName, outFastaName, startIndex=None, stopIndex=None, test=False):
	with open(inFastaName) as filehandle:
		numReads = seqFileToTrimmedFasta_filehandle(filehandle, outFastaName,'fasta',startIndex,stopIndex,test)

	return numReads

def seqFileToTrimmedFasta_filehandle(infile, fastaFileName, seqFileType, startIndex=None, stopIndex=None, test=False):
	curRead = 0
	with open(fastaFileName,'w') as outfile:
		for seq_record in SeqIO.parse(infile,seqFileType):
			outfile.write('>' + seq_record.id + '\n' + str(seq_record.seq)[startIndex:stopIndex] + '\n')
			curRead += 1

			#allow test runs using only the first 100 reads from the fastq file
			if test and curRead >= testLines:
				break

	return curRead

### Bowtie mapping trimmed fasta files ###

def prepFileListsForBowtie(basePath, outfileBaseList, bowtieIndex):
	bowtieMapPath = os.path.join(basePath,'bowtiemaps')
	bowtieUnmapPath = os.path.join(bowtieMapPath,'bowtieUnmapped')
	makeDirectory(bowtieUnmapPath)

	bowtieBaseList = [outfileName + '_' + os.path.split(bowtieIndex)[-1] for outfileName in outfileBaseList]
	bowtieMapPathList = [os.path.join(bowtieMapPath,bowtieMapName) + '.map' for bowtieMapName in bowtieBaseList]
	bowtieUnmapPathList = [os.path.join(bowtieUnmapPath,bowtieMapName) + '.unmap.fa' for bowtieMapName in bowtieBaseList]
	
	alreadyProcessedMapsList = [os.path.join(bowtieMapPath, fileName) for fileName in os.listdir(bowtieMapPath)]

	return bowtieMapPath, bowtieUnmapPath, bowtieBaseList, bowtieMapPathList, bowtieUnmapPathList, alreadyProcessedMapsList

def mapTrimmedFastasToLibrary(fastaFilePathList, bowtieMapPathList, bowtieUnmapPathList, bowtieIndex, numProcessors, mapFilesToSkip=None):
	numProcessors = max(numProcessors,1)

	for fastaFile, mapFile, unmapFile in zip(fastaFilePathList, bowtieMapPathList, bowtieUnmapPathList):
		if mapFilesToSkip == None or mapFile not in mapFilesToSkip:
			printNow('Starting ' + mapFile + ' ...')
			subprocess.call('bowtie -v 0 -a -p ' + repr(numProcessors) + ' --norc ' + bowtieIndex + ' --suppress 2,4,5,6,7,8 -f ' + fastaFile + ' ' + mapFile + ' --un ' + unmapFile, shell=True)


### Map File to Counts File Functions ###

#adapted from Martin's align_primary.py
def countReadsInMapFile(mapFileName, countsFileName, ambigFileName):
	mapCounts = dict()
	ambigReads = set()
	with open(mapFileName) as infile:
		infilereader = csv.reader(infile, delimiter='\t')
		prevRead = ''
		for line in infilereader:
			readId, mapId = line
			
			if mapId not in mapCounts:
				mapCounts[mapId] = 0

			mapCounts[mapId] += 1

			if readId == prevRead:
				ambigReads.add(readId)
			prevRead = readId

	with open(countsFileName,'w') as outfile:
		outfilewriter = csv.writer(outfile, delimiter = '\t')
		for mapId in sorted(mapCounts.keys()):
			outfilewriter.writerow([mapId, mapCounts[mapId]])

	with open(ambigFileName,'w') as outfile:
		for read in sorted(list(ambigReads)):
			outfile.write(read + '\n')

	return countsFileName, sum(mapCounts.values()), len(ambigReads)

def countReadsInMapFileWrapper(args):
	return countReadsInMapFile(*args)


### Utility Functions ###
def parseSeqFileNames(fileNameList):
	infileList = []
	outfileBaseList = []

	for inputFileName in fileNameList:					#iterate through entered filenames for sequence files
		for filename in glob.glob(inputFileName): 		#generate all possible files given wildcards
			for fileType in zip(*acceptedFileTypes)[0]:	#iterate through allowed filetypes
				if fnmatch.fnmatch(filename,fileType):
					infileList.append(filename)
					outfileBaseList.append(os.path.split(filename)[-1].split('.')[0])

	return infileList, outfileBaseList

def makeDirectory(path):
	try:
		os.makedirs(path)
	except OSError:
		#printNow(path + ' already exists')
		pass

def printNow(printInput):
	print printInput
	sys.stdout.flush()

### Global variables ###
acceptedFileTypes = [('*.fastq.gz',fastqGzToTrimmedFasta),
					('*.fastq', fastqToTrimmedFasta),
					('*.fa', fastaToTrimmedFasta),
					('*.fasta', fastaToTrimmedFasta),
					('*.fna', fastaToTrimmedFasta)]


testLines = 100


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process raw sequencing data from screens to counts files in parallel')
	parser.add_argument('Bowtie_Index', help='Bowtie index root name (without .1.ebwt).')
	parser.add_argument('Out_File_Path')
	parser.add_argument('Seq_File_Names', nargs='+', help='Name(s) of sequencing file(s). Unix wildcards can be used to select multiple files at once. The script will search for all *.fastq.gz, *.fastq, and *.fa(/fasta/fna) files with the given wildcard name.')
			
	parser.add_argument('-u','--Unaligned_Indices', nargs='+', help='Bowtie indices to test unaligned reads for possible cross-contaminations.')
	parser.add_argument('-p','--processors', type=int, default = 1)
	parser.add_argument('--trim_start', type=int)
	parser.add_argument('--trim_end', type=int)
	parser.add_argument('--test', action='store_true', default=False, help='Run the entire script on only the first %d reads of each file. Be sure to delete or move all test files before re-running script as they will not be overwritten.' % testLines)

	args = parser.parse_args()
	#printNow(args)

	numProcessors = max(args.processors, 1)

	infileList, outfileBaseList = parseSeqFileNames(args.Seq_File_Names)
			
	## unzip and trim all files ##
	trimmedFastaPath = os.path.join(args.Out_File_Path,'trimmed_fastas')
	makeDirectory(trimmedFastaPath)

	fastaFileNameList = [outfileName + '_trim.fa' for outfileName in outfileBaseList] #now that fasta files are allowed, add _trim to avoid clash if seq file path and out path are same
	fastaFilePathList = [os.path.join(trimmedFastaPath, fastaFileName) for fastaFileName in fastaFileNameList]
	allFilesToProcess = [fastaFileName for fastaFileName in zip(infileList,fastaFileNameList) if fastaFileName[1] not in os.listdir(trimmedFastaPath)]
	if len(allFilesToProcess) != 0:
		infilesToProcess, fastasToProcess = zip(*allFilesToProcess)
	else:
		infilesToProcess, fastasToProcess = [],[]
	fastaFilePathsToProcess = [os.path.join(trimmedFastaPath, fastaFileName) for fastaFileName in fastasToProcess]

	
	if len(fastaFilePathsToProcess) != 0:
		printNow('Parsing and trimming sequence files...')
		sys.stdout.flush()

		pool = multiprocessing.Pool(min(len(infilesToProcess),numProcessors))

		resultList = parallelSeqFileToTrimmedFasta(infilesToProcess, fastaFilePathsToProcess,args.Out_File_Path,pool, args.trim_start, args.trim_end, args.test)
		for result in resultList:
			print result[0] + ':\t' + repr(result[1]) + ' reads'
		
		pool.close()
		pool.join()
		
	printNow('Done parsing and trimming sequence files')

	## bowtie map all files ##
	bowtieIndex = args.Bowtie_Index
	printNow('Mapping all sequencing runs to index: ' + bowtieIndex)

	bowtieMapPath, bowtieUnmapPath, bowtieBaseList, bowtieMapPathList, bowtieUnmapPathList, alreadyProcessedMapsList = \
		prepFileListsForBowtie(args.Out_File_Path, outfileBaseList, bowtieIndex)
	
	mapTrimmedFastasToLibrary(fastaFilePathList, bowtieMapPathList, bowtieUnmapPathList, bowtieIndex, numProcessors, alreadyProcessedMapsList)

	if args.Unaligned_Indices != None:
		prevUnmapPath = bowtieUnmapPath
		prevbowtieUnmapPathList = bowtieUnmapPathList
		for unalignedIndex in args.Unaligned_Indices:
			printNow('Mapping remaining unaligned reads to index: ' + bowtieIndex)

			curbowtieMapPath, curbowtieUnmapPath, curbowtieBaseList, curbowtieMapPathList, curbowtieUnmapPathList, curalreadyProcessedMapsList = \
				prepFileListsForBowtie(prevUnmapPath, [os.path.split(unmapPathName)[-1] for unmapPathName in prevbowtieUnmapPathList], unalignedIndex)

			mapTrimmedFastasToLibrary(prevbowtieUnmapPathList, curbowtieMapPathList, curbowtieUnmapPathList, unalignedIndex, numProcessors, curalreadyProcessedMapsList)

			prevUnmapPath = curbowtieUnmapPath
			prevbowtieUnmapPathList = curbowtieUnmapPathList


	printNow('Done bowtie mapping')

	## count bowtie read mappings ##
	countFilePath = os.path.join(args.Out_File_Path,'count_files')
	ambiguousFilePath = os.path.join(countFilePath, 'ambiguous_reads')
	makeDirectory(ambiguousFilePath)

	countPathList = [os.path.join(countFilePath,bowtieBase + '.counts') for bowtieBase in bowtieBaseList]
	ambigPathList = [os.path.join(ambiguousFilePath,bowtieBase + '.ambig') for bowtieBase in bowtieBaseList]

	filesToProcess = [filePaths for filePaths in zip(bowtieMapPathList,countPathList,ambigPathList) if os.path.split(filePaths[1])[-1] not in os.listdir(countFilePath)]
	#printNow(filesToProcess)

	if len(filesToProcess) != 0:
		printNow('Counting mapped reads...')

		pool = multiprocessing.Pool(min(len(filesToProcess), numProcessors))
		
		countResults = pool.map(countReadsInMapFileWrapper, filesToProcess)
		for result in countResults:
			printNow('%s: \t%e total mappings\t%e ambiguous reads' % result)

		pool.close()
		pool.join()
		
	printNow('Done counting mapped reads')