# Planned functionality:
# merge counts files into a data table, combine reads from multiple sequencing runs,
#  filter by read counts, generate phenotype scores, average replicates

import pandas as pd
import os
import sys
from Bio import Seq, Alphabet
import numpy as np
import scipy as sp
import fnmatch
import scipy.stats.mstats as ms
import matplotlib.pyplot as plt

from expt_config_parser import parseExptConfig, parseLibraryConfig

defaultLibConfigName = 'library_config.txt'

#a screen processing pipeline that requires just a config file and a directory of supported libraries
#error checking in config parser is fairly robust, so not checking for input errors here
def processExperimentsFromConfig(configFile, libraryDirectory):
    #load in the supported libraries and sublibraries
    librariesToSublibraries, librariesToTables = parseLibraryConfig(os.path.join(libraryDirectory, defaultLibConfigName))

    exptParameters, parseStatus, parseString = parseExptConfig(configFile, librariesToSublibraries)

    print parseString
    sys.stdout.flush()

    if parseStatus > 0: #Critical errors in parsing
        print 'Exiting due to parsing errors\n'
        return

    outbase = os.path.join(exptParameters['output_folder'],exptParameters['experiment_name'])

    #load in library table and filter to requested sublibraries
    print 'Accessing library information'
    sys.stdout.flush()

    libraryTable = pd.read_csv(os.path.join(libraryDirectory, librariesToTables[exptParameters['library']]), sep = '\t', tupleize_cols=False, header=0, index_col=0).sort_index()
    sublibColumn = libraryTable.apply(lambda row: row['sublibrary'] in exptParameters['sublibraries'], axis=1)

    libraryTable[sublibColumn].to_csv(outbase + '_librarytable.txt', sep='\t', tupelize_cols = False)

    #load in counts, create table of total counts in each and each file as a column
    print 'Loading counts data'
    sys.stdout.flush()

    columnDict = dict()
    for tup in sorted(exptParameters['counts_file_list']):
        if tup in columnDict:
            print 'Asserting that tuples of condition, replicate, and count file should be unique; are the cases where this should not be enforced?'
            raise Exception('condition, replicate, and count file combination already assigned')
        
        countSeries = readCountsFile(tup[2]).reset_index().drop_duplicates('id').set_index('id') #for now also dropping duplicate ids in counts for overlapping linc sublibraries
        countSeries = libraryTable[sublibColumn].align(countSeries, axis=0, join='left', fill_value=0)[1] #expand series to fill 0 for every missing entry

        columnDict[tup] = countSeries['counts'] #[sublibColumn] #then shrink series to only desired sublibraries

    # print columnDict
    countsTable = pd.DataFrame(columnDict)#, index=libraryTable[sublibColumn].index)
    countsTable.to_csv(outbase + '_rawcountstable.txt', sep='\t', tupleize_cols = False)
    countsTable.sum().to_csv(outbase + '_rawcountstable_summary.txt', sep='\t')

    #merge counts for same conditions/replicates, and create summary table
    #save scatter plot before each merger, and histogram of counts post mergers
    print 'Merging experiment counts split across lanes/indexes'
    print 'Generating scatter plots of counts pre-merger and histograms of counts post-merger'
    sys.stdout.flush()

    exptGroups = countsTable.groupby(level=[0,1], axis=1)
    for (condition, replicate), countsCols in exptGroups:
        if len(countsCols.columns) == 1:
            continue

        for i, col1 in enumerate(countsCols):
            for j, col2 in enumerate(countsCols):
                if j > i: #enforce that each pair is compared once
                    rasteredScatter(countsCols[col1],countsCols[col2],'\n'.join(col1),
                        '\n'.join(col2),outbase + '_premergescatter_%s_%s_%dv%d.svg' % (condition,replicate,i,j))

    mergedCountsTable = exptGroups.aggregate(np.sum)
    mergedCountsTable.to_csv(outbase + '_mergedcountstable.txt', sep='\t', tupleize_cols = False)
    mergedCountsTable.sum().to_csv(outbase + '_mergedcountstable_summary.txt', sep='\t')

    for col in mergedCountsTable:
        generateHistogram(mergedCountsTable[col],', '.join(col),outbase + '_postmergehist_%s_%s.svg' % (condition,replicate))

    #create pairs of columns for each comparison, filter to na, then generate sgRNA phenotype score
    print 'Computing sgRNA phenotype scores'
    sys.stdout.flush()

    growthValueDict = {(tup[0],tup[1]):tup[2] for tup in exptParameters['growth_value_tuples']}
    phenotypeList = list(set(zip(*exptParameters['condition_tuples'])[0]))
    replicateList = list(set(zip(*exptParameters['counts_file_list'])[1]))

    phenotypeScoreDict = dict()
    for (phenotype, condition1, condition2) in exptParameters['condition_tuples']:
        for replicate in replicateList:
            column1 = mergedCountsTable[(condition1,replicate)]
            column2 = mergedCountsTable[(condition2,replicate)]
            filtCols = filterLowCounts(pd.concat((column1, column2), axis = 1), exptParameters['filter_type'], exptParameters['minimum_reads'])
            

            score = computePhenotypeScore(filtCols[(condition1, replicate)], filtCols[(condition2,replicate)], 
                libraryTable[sublibColumn], growthValueDict[(phenotype,replicate)], 
                exptParameters['pseudocount_behavior'], exptParameters['pseudocount'])

            phenotypeScoreDict[(phenotype,replicate)] = score

    
    #scatterplot sgRNAs for all replicates, then average together and add columns to phenotype score table
    if len(replicateList) > 1:
        print 'Plotting and averaging replicates'
        sys.stdout.flush()

        for phenotype in phenotypeList:
            for i, rep1 in enumerate(replicateList):
                for j, rep2 in enumerate(replicateList):
                    if j > i:
                        rasteredScatter(phenotypeScoreDict[(phenotype,rep1)],phenotypeScoreDict[(phenotype,rep2)],
                            ', '.join((phenotype,rep1)), ', '.join((phenotype,rep2)),
                            outbase + '_phenotypescatter_%s_%sv%s.svg' % (condition,rep1,rep2))

            repCols = pd.DataFrame({(phen,rep):col for (phen,rep), col in phenotypeScoreDict.iteritems() if phen == phenotype})
            phenotypeScoreDict[(phenotype,'ave_' + '_'.join(replicateList))] = repCols.mean(axis=1,skipna=False) #average nan and real to nan; otherwise this could lead to data points with just one rep informing results

    phenotypeTable = pd.DataFrame(phenotypeScoreDict)
    phenotypeTable.to_csv(outbase + '_phenotypetable.txt', sep='\t', tupleize_cols = False)

    #generate pseudogenes
    negTable = phenotypeTable.loc[libraryTable[sublibColumn].loc[:,'gene'] == 'negative_control',:]

    if exptParameters['generate_pseudogene_dist'] == True and len(exptParameters['analyses']) > 0:
        print 'Generating a pseudogene distribution from negative controls'
        sys.stdout.flush()

        pseudoTableList = []
        pseudoLibTables = []
        negValues = negTable.values
        negColumns = negTable.columns
        for pseudogene in range(exptParameters['num_pseudogenes']):
            randIndices = np.random.randint(0, len(negTable), exptParameters['pseudogene_size'])
            pseudoTable = negValues[randIndices,:]
            pseudoIndex = ['pseudo_%d_%d' % (pseudogene,i) for i in range(exptParameters['pseudogene_size'])]
            pseudoSeqs = ['seq_%d_%d' % (pseudogene,i) for i in range(exptParameters['pseudogene_size'])] #so pseudogenes aren't treated as duplicates
            pseudoTableList.append(pd.DataFrame(pseudoTable,index=pseudoIndex,columns=negColumns))
            pseudoLib = pd.DataFrame({'gene':['pseudo_%d'%pseudogene]*exptParameters['pseudogene_size'],
                'transcripts':['na']*exptParameters['pseudogene_size'],
                'sequence':pseudoSeqs},index=pseudoIndex)
            pseudoLibTables.append(pseudoLib)

        phenotypeTable = phenotypeTable.append(pd.concat(pseudoTableList))
        libraryTableGeneAnalysis = libraryTable[sublibColumn].append(pd.concat(pseudoLibTables))
    else:
        libraryTableGeneAnalysis = libraryTable[sublibColumn]

    #return phenotypeTable, libraryTable

    #compute gene scores for replicates, averaged reps, and pseudogenes
    if len(exptParameters['analyses']) > 0:
        print 'Computing gene scores'
        sys.stdout.flush()

        phenotypeTable_deduplicated = phenotypeTable.loc[libraryTableGeneAnalysis.drop_duplicates(['gene','sequence']).index]
        if exptParameters['collapse_to_transcripts'] == True:
            geneGroups = phenotypeTable_deduplicated.loc[libraryTableGeneAnalysis.loc[:,'gene'] != 'negative_control',:].groupby([libraryTableGeneAnalysis['gene'],libraryTableGeneAnalysis['transcripts']])
        else:
            geneGroups = phenotypeTable_deduplicated.loc[libraryTableGeneAnalysis.loc[:,'gene'] != 'negative_control',:].groupby(libraryTableGeneAnalysis['gene'])

        analysisTables = []
        for analysis in exptParameters['analyses']:
            print '--' + analysis
            sys.stdout.flush()

            analysisTables.append(applyGeneScoreFunction(geneGroups, negTable, analysis, exptParameters['analyses'][analysis]))

        geneTable = pd.concat(analysisTables, axis=1).reorder_levels([1,2,0],axis=1).sort_index(axis=1)
        geneTable.to_csv(outbase + '_genetable.txt',sep='\t', tupleize_cols = False)

        ### collapse the gene-transcript indices into a single score for a gene by best MW p-value, where applicable
        if exptParameters['collapse_to_transcripts'] == True and 'calculate_mw' in exptParameters['analyses']:
            print 'Collapsing transcript scores to gene scores'
            sys.stdout.flush()

            geneTableCollapsed = scoreGeneByBestTranscript(geneTable)
            geneTableCollapsed.to_csv(outbase + '_genetable_collapsed.txt',sep='\t', tupleize_cols = False)


    #generate summary graphs depending on which analyses were selected


    print 'Done!'

#given a gene table indexed by both gene and transcript, score genes by the best m-w p-value per phenotype/replicate
def scoreGeneByBestTranscript(geneTable):
    geneTableTransGroups = geneTable.reorder_levels([2,0,1],axis=1)['Mann-Whitney p-value'].reset_index().groupby('gene')

    bestTranscriptFrame = geneTableTransGroups.apply(getBestTranscript)

    tupList = []
    bestTransList = []
    for tup, group in geneTable.groupby(level=range(2),axis=1):
        tupList.append(tup)
        curFrame = geneTable.loc[zip(bestTranscriptFrame.index,bestTranscriptFrame[tup]),tup]
        bestTransList.append(curFrame.reset_index().set_index('gene'))

    return pd.concat(bestTransList, axis=1, keys=tupList)

def getBestTranscript(group):
    #set the index to be transcripts and then get the index with the lowest p-value for each cell
    return group.set_index('transcripts').drop(('gene',''),axis=1).idxmin() 

#an example processing pipeline, standardized for crispricin screens
def processCrispricinExperiments(outbase, libraryFastaFileName,experimentFileName, gkFileName, filterThreshold = 25, filterSet = None, normToNegs = True):
    libraryTable = readLibraryFile(libraryFastaFileName,crisprElementType,crisprGeneName,[crisprGuideLength_v1, crisprGuideSequence_rctrimmed28])
    print 'Merging counts for each experiment...'
    rawCounts, allExpts = mergeCountsForExperiments(experimentFileName,libraryTable)
    print 'Filtering guides by low read counts across all conditions/replicates...\n # guides filtered:'
    filteredCounts = filterCountsPerExperiment(filterThreshold, allExpts, libraryTable).sort(axis=1)
    if filterSet != None:
        filteredCounts = filteredCounts.loc[set(filteredCounts.index.values) - filterSet].align(filteredCounts,axis=0)[0]
    
    #table1 = pd.concat([libraryTable, rawCounts, allExpts, filteredCounts], axis=1, keys = ['library_table', 'raw_counts', 'expt_counts', 'filtered_counts'])
    
    libraryTable.to_csv(outbase+'_librarytable.txt', sep='\t', tupelize_cols = False)
    rawCounts.to_csv(outbase+'_rawcounts.txt', sep='\t', tupelize_cols = False)
    allExpts.sort(axis=1).to_csv(outbase+'_mergedexptcounts.txt', sep='\t', tupelize_cols = False)
    filteredCounts.to_csv(outbase+'_filteredcounts.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_countstable.txt',sep='\t', tupleize_cols=False, header=range(3), index_col=0)

    print 'Computing phenotype scores...'
    gkdict = parseGKFile(gkFileName)
    gammas, taus, rhos = computeAllPhenotypeScores('t0','cycled','ricin', filteredCounts, libraryTable, gkdict, normToNegs=normToNegs)
    mergedScores = pd.concat([gammas,taus,rhos], axis=1).sort(axis=1)
    aveScores = averagePhenotypeScores(mergedScores).sort(axis=1)
    
    #table2 = pd.concat([libraryTable, mergedScores,aveScores], axis = 1, keys = ['library_table','replicate_scores', 'averaged_scores'])
    mergedScores.to_csv(outbase+'_replicatephenotypescores.txt', sep='\t', tupelize_cols = False)
    aveScores.to_csv(outbase+'_averagedphenotypescores.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_phenotypescores.txt',sep='\t', tupleize_cols=False, header=range(3), index_col=0)

    print 'Computing gene scores...'
    geneScores = computeGeneScores(libraryTable, mergedScores, normToNegs=normToNegs).sort(axis=1)
    aveGeneScores = computeGeneScores(libraryTable, aveScores, normToNegs=normToNegs).sort(axis=1)
    
    #table3 = pd.concat([geneScores, aveGeneScores], axis = 1, keys = ['replicate_gene_pvals','averaged_gene_pvals'])
    geneScores.to_csv(outbase+'_replicategenescores.txt', sep='\t', tupelize_cols = False)
    aveGeneScores.to_csv(outbase+'_averagedgenescores.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_genescores.txt',sep='\t', tupleize_cols=False, header=range(4), index_col=0)
    print 'Done!'

#an example processing pipeline, standardized for crispricin screens
def processEssentialExperiments(outbase, libraryFastaFileName,experimentFileName, gkFileName, filterThreshold = 25, filterSet = None, normToNegs = True):
    libraryTable = readLibraryFile(libraryFastaFileName,crisprElementTypeSublibraries,crisprGeneNameSublibraries,[])#[crisprGuideLength_v2, crisprGuideSequence_trimmed35])
    
    print 'Merging counts for each experiment...'
    rawCounts, allExpts = mergeCountsForExperiments(experimentFileName,libraryTable)
    print 'Filtering guides by low read counts across all conditions/replicates...\n # guides filtered:'
    filteredCounts = filterCountsPerExperiment(filterThreshold, allExpts, libraryTable).sort(axis=1)
    if filterSet != None:
        filteredCounts = filteredCounts.loc[set(filteredCounts.index.values) - filterSet].align(filteredCounts,axis=0)[0]
    
    #table1 = pd.concat([libraryTable, rawCounts, allExpts, filteredCounts], axis=1, keys = ['library_table', 'raw_counts', 'expt_counts', 'filtered_counts'])
    
    libraryTable.to_csv(outbase+'_librarytable.txt', sep='\t', tupelize_cols = False)
    rawCounts.to_csv(outbase+'_rawcounts.txt', sep='\t', tupelize_cols = False)
    allExpts.sort(axis=1).to_csv(outbase+'_mergedexptcounts.txt', sep='\t', tupelize_cols = False)
    filteredCounts.to_csv(outbase+'_filteredcounts.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_countstable.txt',sep='\t', tupleize_cols=False, header=range(3), index_col=0)

    print 'Computing phenotype scores...'
    gkdict = parseGKFile(gkFileName)
    gammas, taus, rhos = computeAllPhenotypeScores('T0','cycled','rigo', filteredCounts, libraryTable, gkdict, normToNegs=normToNegs)
    mergedScores = pd.concat([gammas,taus,rhos], axis=1).sort(axis=1)
    aveScores = averagePhenotypeScores(mergedScores).sort(axis=1)
    
    #table2 = pd.concat([libraryTable, mergedScores,aveScores], axis = 1, keys = ['library_table','replicate_scores', 'averaged_scores'])
    mergedScores.to_csv(outbase+'_replicatephenotypescores.txt', sep='\t', tupelize_cols = False)
    aveScores.to_csv(outbase+'_averagedphenotypescores.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_phenotypescores.txt',sep='\t', tupleize_cols=False, header=range(3), index_col=0)

    #print 'Computing gene scores...'
    #geneScores = computeGeneScores(libraryTable, mergedScores, normToNegs=normToNegs).sort(axis=1)
    #aveGeneScores = computeGeneScores(libraryTable, aveScores, normToNegs=normToNegs).sort(axis=1)
    
    #table3 = pd.concat([geneScores, aveGeneScores], axis = 1, keys = ['replicate_gene_pvals','averaged_gene_pvals'])
    #geneScores.to_csv(outbase+'_replicategenescores.txt', sep='\t', tupelize_cols = False)
    #aveGeneScores.to_csv(outbase+'_averagedgenescores.txt', sep='\t', tupelize_cols = False)
    #testReadback = pd.read_csv('testfiles/testtable_genescores.txt',sep='\t', tupleize_cols=False, header=range(4), index_col=0)
    print 'Done!'

#return Series of counts from a counts file indexed by element id
def readCountsFile(countsFileName):
    countsTable = pd.read_csv(countsFileName, header=None, delimiter='\t', names=['id','counts'])
    countsTable.index = countsTable['id']
    return countsTable['counts']


#return DataFrame of library features indexed by element id
def readLibraryFile(libraryFastaFileName, elementTypeFunc, geneNameFunc, miscFuncList=None):
    elementList = []
    with open(libraryFastaFileName) as infile:
        idLine = infile.readline()
        while idLine != '':
            seqLine = infile.readline()
            if idLine[0] != '>' or seqLine == None:
                raise ValueError('Error parsing fasta file')
                
            elementList.append((idLine[1:].strip(), seqLine.strip()))
            
            idLine = infile.readline()

    elementIds, elementSeqs = zip(*elementList)
    libraryTable = pd.DataFrame(np.array(elementSeqs), index=np.array(elementIds), columns=['aligned_seq'], dtype='object')

    libraryTable['element_type'] = elementTypeFunc(libraryTable)
    libraryTable['gene_name'] = geneNameFunc(libraryTable)
    
    if miscFuncList != None:
        colList = [libraryTable]
        for miscFunc in miscFuncList:
            colList.append(miscFunc(libraryTable))
        if len(colList) != 1:
            libraryTable = pd.concat(colList, axis=1)

    return libraryTable

#print all counts file paths, to assist with making an experiment table
def printCountsFilePaths(baseDirectoryPathList):
    print 'Make a tab-delimited file with the following columns:'
    print 'counts_file\texperiment\tcondition\treplicate_id'
    print 'and the following list in the counts_file column:'
    for basePath in baseDirectoryPathList:
        for root, dirs, filenames in os.walk(basePath):
            for filename in fnmatch.filter(filenames,'*.counts'):
                print os.path.join(root, filename)

def mergeCountsForExperiments(experimentFileName, libraryTable):
    exptTable = pd.read_csv(experimentFileName, delimiter='\t')
    print exptTable

    # load in all counts independently
    countsCols = []
    for countsFile in exptTable['counts_file']:
        countsCols.append(readCountsFile(countsFile))

    countsTable = pd.concat(countsCols, axis=1, keys=exptTable['counts_file']).align(libraryTable,axis=0)[0]
    
    countsTable = countsTable.fillna(value = 0) #nan values are 0 values, will use nan to filter out elements later

    #print countsTable.head()
    
    # convert counts columns to experiments, summing when reads across multiple lanes
    exptTuples = [(exptTable.loc[row,'experiment'],exptTable.loc[row,'condition'],exptTable.loc[row,'replicate_id']) for row in exptTable.index]
    exptTuplesToRuns = dict()
    for i, tup in enumerate(exptTuples):
        if tup not in exptTuplesToRuns:
            exptTuplesToRuns[tup] = []
        exptTuplesToRuns[tup].append(exptTable.loc[i,'counts_file'])

    #print exptTuplesToRuns

    exptColumns = []
    for tup in sorted(exptTuplesToRuns.keys()):
        if len(exptTuplesToRuns[tup]) == 1:
            exptColumns.append(countsTable[exptTuplesToRuns[tup][0]])
        else:
            column = countsTable[exptTuplesToRuns[tup][0]]
            for i in range(1,len(exptTuplesToRuns[tup])):
                column += countsTable[exptTuplesToRuns[tup][i]]

            exptColumns.append(column)

    #print len(exptColumns), exptColumns[-1]

    exptsTable = pd.concat(exptColumns, axis = 1, keys=sorted(exptTuplesToRuns.keys()))
    exptsTable.columns = pd.MultiIndex.from_tuples(sorted(exptTuplesToRuns.keys()))
    #print exptsTable

    #mergedTable = pd.concat([libraryTable,countsTable,exptsTable],axis=1, keys = ['library_properties','raw_counts', 'merged_experiments'])

    return countsTable, exptsTable

#filter out reads if /all/ reads for an expt accross replicates/conditions < min_reads
def filterCountsPerExperiment(min_reads, exptsTable,libraryTable):
    experimentGroups = []
    
    exptTuples = exptsTable.columns

    exptSet = set([tup[0] for tup in exptTuples])
    for expt in exptSet:
        exptDf = exptsTable[[tup for tup in exptTuples if tup[0] == expt]]
        exptDfUnderMin = (exptDf < min_reads).all(axis=1)
        exptDfFiltered = exptDf.align(exptDfUnderMin[exptDfUnderMin == False], axis=0, join='right')[0]
        experimentGroups.append(exptDfFiltered)
        
        print expt, len(exptDfUnderMin[exptDfUnderMin == True])

    resultTable = pd.concat(experimentGroups, axis = 1).align(libraryTable, axis=0)[0]

    return resultTable

#more flexible read filtering 
#keep row if either both/all columns are above threshold, or if either/any column is
#in other words, mask if any column is below threshold or only if all columns are below
def filterLowCounts(countsColumns, filterType, filterThreshold):
    if filterType == 'both' or filterType == 'all':
        failFilterColumn = countsColumns.apply(lambda row: min(row) < filterThreshold, axis = 1)
    elif filterType == 'either' or filterType == 'any':
        failFilterColumn = countsColumns.apply(lambda row: max(row) < filterThreshold, axis = 1)
    else:
        raise ValueError('filter type not recognized or not implemented')

    resultTable = countsColumns.copy()
    resultTable.loc[failFilterColumn,:] = np.nan

    return resultTable


#compute phenotype scores for each experiment
#default pseudocount behavior is +1 for any element with a zero value
def computeAllPhenotypeScores(startCondition, endUntreatedCondition, endTreatedCondition, exptsTable, libraryTable, exptToGKvalues, pseudocounts = 'default', normToNegs=True):
    gammaList = []
    tauList = []
    rhoList = []

    exptTuples = exptsTable.columns
    exptsToReplicates = dict()
    for tup in exptTuples:
        if tup[0] not in exptsToReplicates:
            exptsToReplicates[tup[0]] = set()
        exptsToReplicates[tup[0]].add(tup[2])

    labels = []
    for expt in exptsToReplicates:
        for rep in exptsToReplicates[expt]:
            labels.append((expt,rep))
            gammaList.append(computePhenotypeScore(exptsTable[(expt,startCondition,rep)], \
                exptsTable[(expt,endUntreatedCondition,rep)], libraryTable, exptToGKvalues[(expt,rep)][0],pseudocounts, normToNegs))
            tauList.append(computePhenotypeScore(exptsTable[(expt,startCondition,rep)], \
                exptsTable[(expt,endTreatedCondition,rep)], libraryTable, exptToGKvalues[(expt,rep)][0] - exptToGKvalues[(expt,rep)][1],pseudocounts, normToNegs))
            rhoList.append(computePhenotypeScore(exptsTable[(expt,endUntreatedCondition,rep)], \
                exptsTable[(expt,endTreatedCondition,rep)], libraryTable, exptToGKvalues[(expt,rep)][1],pseudocounts, normToNegs))

    #print labels
    gammas = pd.concat(gammaList, axis = 1, keys = [(lab[0],'gamma',lab[1]) for lab in labels]).align(libraryTable, axis = 0)[0]
    gammas.columns = pd.MultiIndex.from_tuples([(lab[0],'gamma',lab[1]) for lab in labels])
    taus = pd.concat(tauList, axis = 1, keys = [(lab[0],'tau',lab[1]) for lab in labels]).align(libraryTable, axis = 0)[0]
    taus.columns = pd.MultiIndex.from_tuples([(lab[0],'tau',lab[1]) for lab in labels])
    rhos = pd.concat(rhoList, axis = 1, keys = [(lab[0],'rho',lab[1]) for lab in labels]).align(libraryTable, axis = 0)[0]
    rhos.columns = pd.MultiIndex.from_tuples([(lab[0],'rho',lab[1]) for lab in labels])

    return gammas, taus, rhos

#compute phenotype scores for any given comparison of two conditions
#edited 11/13/2014, so computeAllPhenotypeScores and example pipelines may require corrections
def computePhenotypeScore(counts1, counts2, libraryTable, growthValue, pseudocountBehavior, pseudocountValue, normToNegs=True):
    combinedCounts = pd.concat([counts1,counts2],axis = 1)

    #pseudocount
    if pseudocountBehavior == 'default' or pseudocountBehavior == 'zeros only':
        defaultBehavior = lambda row: row if min(row) != 0 else row + pseudocountValue
        combinedCountsPseudo = combinedCounts.apply(defaultBehavior, axis = 1)
    elif pseudocountBehavior == 'all values':
        combinedCountsPseudo = combinedCounts.apply(lambda row: row + pseudocountValue, axis = 1)
    elif pseudocountBehavior == 'filter out':
        combinedCountsPseudo = combinedCounts.copy()
        zeroRows = combinedCounts.apply(lambda row: min(row) <= 0, axis = 1)
        combinedCountsPseudo.loc[zeroRows,:] = np.nan
    else:
        raise ValueError('Pseudocount behavior not recognized or not implemented')

    totalCounts = combinedCountsPseudo.sum()
    countsRatio = float(totalCounts[0])/totalCounts[1]

    #get neg control log2e--does this need WTlog2E??
    if normToNegs == True:
        negCounts = combinedCountsPseudo.align(libraryTable[libraryTable['gene'] == 'negative_control'],axis=0,join='inner')[0]
        #print negCounts
    else:
        negCounts = combinedCountsPseudo
    neglog2e = negCounts.apply(calcLog2e, countsRatio=countsRatio, growthValue=1, wtLog2E=0, axis=1).median() #no growth value used in martin's
    #print neglog2e
    #compute scores
    scores = combinedCountsPseudo.apply(calcLog2e, countsRatio=countsRatio, growthValue=growthValue, wtLog2E=neglog2e, axis=1)

    return scores

def calcLog2e(row, countsRatio, growthValue, wtLog2E):
    return (np.log2(countsRatio*row[1]/row[0]) - wtLog2E) / growthValue

#average replicate phenotype scores
def averagePhenotypeScores(scoreTable):

    exptTuples = scoreTable.columns
    exptsToReplicates = dict()
    for tup in exptTuples:
        if (tup[0],tup[1]) not in exptsToReplicates:
            exptsToReplicates[(tup[0],tup[1])] = set()
        exptsToReplicates[(tup[0],tup[1])].add(tup[2])

    averagedColumns = []
    labels = []
    for expt in exptsToReplicates:
        exptDf = scoreTable[[(expt[0],expt[1],rep_id) for rep_id in exptsToReplicates[expt]]]
        averagedColumns.append(exptDf.mean(axis=1))
        labels.append((expt[0],expt[1],'ave_'+'_'.join(exptsToReplicates[expt])))

    resultTable = pd.concat(averagedColumns, axis = 1, keys=labels).align(scoreTable, axis=0)[0]
    resultTable.columns = pd.MultiIndex.from_tuples(labels)

    return resultTable

def computeGeneScores(libraryTable, scoreTable, normToNegs = True):
    geneGroups = scoreTable.groupby(libraryTable['gene_name'])
    
    scoredColumns = []
    for expt in scoreTable.columns:
        if normToNegs == True:
            negArray = np.ma.array(data=scoreTable[expt].loc[geneGroups.groups['negative_control']].dropna(),mask=False)
        else:
            negArray = np.ma.array(data=scoreTable[expt].dropna(),mask=False)
        
        colList = []
        groupList = []
        for name, group in geneGroups:
            if name == 'negative_control':
                continue
            colList.append(geneStats(group[expt],negArray)) #group[expt].apply(geneStats, axis = 0, negArray = negArray))
            groupList.append(name)
            
        scoredColumns.append(pd.DataFrame(np.array(colList), index = groupList, columns = [('KS'),('KS_sign'),('MW')]))
    
    #return scoredColumns
    return pd.concat(scoredColumns, axis = 1, keys=scoreTable.columns)

def geneStats(scoreColumn, negArray):
    scoreArray = np.ma.array(data=scoreColumn.dropna(), mask=False)
    ksPval = ms.ks_twosamp(scoreArray, negArray)[1]
    
    ksHi = ms.ks_twosamp(scoreArray, negArray, alternative = 'less')[1]
    ksLo = ms.ks_twosamp(scoreArray, negArray, alternative = 'greater')[1]
    if ksHi < ksLo:
        ksSign = 'P'
    else:
        ksSign = 'S'

    mwPval = ms.mannwhitneyu(scoreArray, negArray)[1]

    return ksPval, ksSign, mwPval

#apply gene scoring functions to pre-grouped tables of phenotypes
def applyGeneScoreFunction(groupedPhenotypeTable, negativeTable, analysis, analysisParamList):
    if analysis == 'calculate_ave':
        numToAverage = analysisParamList[0]
        if numToAverage <= 0:
            means = groupedPhenotypeTable.aggregate(np.mean)
            counts = groupedPhenotypeTable.count()
            result = pd.concat([means,counts],axis=1,keys=['average of all phenotypes','average of all phenotypes_sgRNAcount'])
        else:
            means = groupedPhenotypeTable.aggregate(lambda x: averageBestN(x, numToAverage))
            counts = groupedPhenotypeTable.count()
            result = pd.concat([means,counts],axis=1,keys=['average phenotype of strongest %d'%numToAverage, 'sgRNA count_avg'])
    elif analysis == 'calculate_mw':
        pvals = groupedPhenotypeTable.aggregate(lambda x: applyMW(x, negativeTable))
        counts = groupedPhenotypeTable.count()
        result = pd.concat([pvals,counts],axis=1,keys=['Mann-Whitney p-value','sgRNA count_MW'])
    elif analysis == 'calculate_nth':
        nth = analysisParamList[0]
        pvals = groupedPhenotypeTable.aggregate(lambda x: sorted(x, key=abs, reverse=True)[nth-1] if nth <= len(x) else np.nan)
        counts = groupedPhenotypeTable.count()
        result = pd.concat([pvals,counts],axis=1,keys=['%dth best score' % nth,'sgRNA count_nth best'])
    else:
        raise ValueError('Analysis %s not recognized or not implemented' % analysis)

    return result

def averageBestN(column, numToAverage):
    return np.mean(sorted(column.dropna(),key=abs,reverse=True)[:numToAverage])
def applyMW(column, negativeTable):
    if column.count() == 0:
        return np.nan
    else:
        return sp.stats.mannwhitneyu(column.dropna().values, negativeTable[column.name].dropna().values)[1] * 2 #stats.mannwhitneyu is one-tailed!!

# unneccesary and variable behavior--pandas count() automatic discounts nans
# def countAfterFilter(column):
#     return len(column.dropna())


#parse a tab-delimited file with column headers: experiment, replicate_id, G_value, K_value (calculated with martin's parse_growthdata.py)
def parseGKFile(gkFileName):
    gkdict = dict()
    
    with open(gkFileName,'rU') as infile:
        for line in infile:
            if line.split('\t')[0] == 'experiment':
                continue
            else:
                linesplit = line.strip().split('\t')
                gkdict[(linesplit[0],linesplit[1])] = (float(linesplit[2]),float(linesplit[3]))

    return gkdict

#converter function for crispr element names to element type
def crisprElementType(libTable):
    idArray = libTable.index.values
    typeList = []
    for elementId in idArray:
        if elementId[:3] == 'neg':
            typeList.append('negative_control')
        elif elementId[:3] == 'mis':
            typeList.append('mismatch')
        else:
            typeList.append('sample')

    return pd.DataFrame(np.array(typeList), index=libTable.index)

#converter function for crispr element names to element type
def crisprElementTypeSublibraries(libTable):
    idArray = libTable.index.values
    typeList = []
    for elementId in idArray:
        sgId = parseSgId(elementId)
        if sgId['gene_name'][:3] == 'neg':
            typeList.append('negative_control')
        else:
            typeList.append(sgId['gene_name'])
    
        #sgId = elementId.split('=')[1]
        #if sgId[:3] == 'neg' or sgId[:4] == 'CTRL':
        #   typeList.append('negative_control')
        #elif sgId[:3] == 'mis':
        #   typeList.append('mismatch')
        #else:
        #   typeList.append('sample')

    return pd.DataFrame(np.array(typeList), index=libTable.index)

#converter function for crispr element names to gene name
def crisprGeneName(libTable):
    idArray = libTable.index.values
    nameList = []
    for elementId in idArray:
        if elementId[:3] == 'neg':
            nameList.append('negative_control')
        else:
            nameList.append(elementId.split('_')[0])

    return pd.DataFrame(np.array(nameList), index=libTable.index)

#converter function for crispr element names to gene name
def crisprGeneNameSublibraries(libTable):
    idArray = libTable.index.values
    nameList = []
    for elementId in idArray:
        sgId = parseSgId(elementId)
        if sgId['gene_name'][:3] == 'neg':
            nameList.append('negative_control')
        else:
            nameList.append(sgId['gene_name'])
            
        #sgId = elementId.split('=')[1]
        #if sgId[:3] == 'neg':  # or sgId[:4] == 'CTRL':
        #   nameList.append('negative_control')
        #else:
        #   if sgId[:2] == 'sg':
        #       nameList.append(sgId[2:].split('_')[0])
        #   else:
        #       nameList.append(sgId.split('_')[0])

    return pd.DataFrame(np.array(nameList), index=libTable.index)

def parseSgId(sgId):
    parseDict = dict()
    
    #sublibrary
    if len(sgId.split('=')) == 2:
        parseDict['Sublibrary'] = sgId.split('=')[0]
        remainingId = sgId.split('=')[1]
    else:
        parseDict['Sublibrary'] = None
        remainingId = sgId
        
    #gene name and strand
    underscoreSplit = remainingId.split('_')
    
    for i,item in enumerate(underscoreSplit):
        if item == '+':
            strand = '+'
            geneName = '_'.join(underscoreSplit[:i])
            remainingId = '_'.join(underscoreSplit[i+1:])
            break
        elif item == '-':
            strand = '-'
            geneName = '_'.join(underscoreSplit[:i])
            remainingId = '_'.join(underscoreSplit[i+1:])
            break
        else:
            continue
            
    parseDict['strand'] = strand
    parseDict['gene_name'] = geneName
        
    #position
    dotSplit = remainingId.split('.')
    parseDict['position'] = int(dotSplit[0])
    remainingId = '.'.join(dotSplit[1:])
    
    #length incl pam
    dashSplit = remainingId.split('-')
    parseDict['length'] = int(dashSplit[0])
    remainingId = '-'.join(dashSplit[1:])
    
    #pass score
    tildaSplit = remainingId.split('~')
    parseDict['pass_score'] = tildaSplit[-1]
    remainingId = '~'.join(tildaSplit[:-1]) #should always be length 1 anyway
    
    #transcripts
    parseDict['transcript_list'] = remainingId.split(',')
    
    return parseDict
    
def rasteredScatter(series1,series2,label1,label2,outfilename):
    # print outfilename
    pass

def generateHistogram(series, label, outfilename):
    pass
    
#converter function for crispr element names to guide length
def crisprGuideLength_v1(libTable):
    idArray = libTable.index.values
    return pd.DataFrame(np.array([int(elementId.split('.')[1]) - 3 for elementId in idArray]), index=libTable.index, columns=['guide_length'])

#converter function for rctrimmed28
def crisprGuideSequence_rctrimmed28(libTable):
    idArray = libTable.index.values
    lengthList = [int(elementId.split('.')[1]) - 3 for elementId in idArray]

    return pd.DataFrame(np.array([Seq.Seq(alignedSeq[:lengthList[i]]).reverse_complement().tostring() for i, alignedSeq in enumerate(libTable['aligned_seq'])]), index=libTable.index, columns=['guide_seq'])


