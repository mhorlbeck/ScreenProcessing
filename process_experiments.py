# merge counts files into a data table, combine reads from multiple sequencing runs,
#  filter by read counts, generate phenotype scores, average replicates

import pandas as pd
import os
import sys
import numpy as np
import scipy as sp
from scipy import stats
import fnmatch
import argparse

from expt_config_parser import parseExptConfig, parseLibraryConfig
from fastqgz_to_counts import makeDirectory, printNow
import screen_analysis

defaultLibConfigName = 'library_config.txt'

#a screen processing pipeline that requires just a config file and a directory of supported libraries
#error checking in config parser is fairly robust, so not checking for input errors here
def processExperimentsFromConfig(configFile, libraryDirectory, generatePlots='png'):
    #load in the supported libraries and sublibraries
    try:
        librariesToSublibraries, librariesToTables = parseLibraryConfig(os.path.join(libraryDirectory, defaultLibConfigName))
    except ValueError as err:
        print ' '.join(err.args)
        return

    exptParameters, parseStatus, parseString = parseExptConfig(configFile, librariesToSublibraries)

    printNow(parseString)

    if parseStatus > 0: #Critical errors in parsing
        print 'Exiting due to experiment config file errors\n'
        return

    makeDirectory(exptParameters['output_folder'])
    outbase = os.path.join(exptParameters['output_folder'],exptParameters['experiment_name'])
    
    if generatePlots != 'off':
        plotDirectory = os.path.join(exptParameters['output_folder'],exptParameters['experiment_name'] + '_plots')
        makeDirectory(plotDirectory)
    
        screen_analysis.changeDisplayFigureSettings(newDirectory=plotDirectory, newImageExtension = generatePlots, newPlotWithPylab = False)
    

    #load in library table and filter to requested sublibraries
    printNow('Accessing library information')

    libraryTable = pd.read_csv(os.path.join(libraryDirectory, librariesToTables[exptParameters['library']]), sep = '\t', tupleize_cols=False, header=0, index_col=0).sort_index()
    sublibColumn = libraryTable.apply(lambda row: row['sublibrary'].lower() in exptParameters['sublibraries'], axis=1)

    if sum(sublibColumn) == 0:
        print 'After limiting analysis to specified sublibraries, no elements are left'
        return

    libraryTable[sublibColumn].to_csv(outbase + '_librarytable.txt', sep='\t', tupleize_cols = False)

    #load in counts, create table of total counts in each and each file as a column
    printNow('Loading counts data')

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
    printNow('Merging experiment counts split across lanes/indexes')
    
    exptGroups = countsTable.groupby(level=[0,1], axis=1)
    mergedCountsTable = exptGroups.aggregate(np.sum)
    mergedCountsTable.to_csv(outbase + '_mergedcountstable.txt', sep='\t', tupleize_cols = False)
    mergedCountsTable.sum().to_csv(outbase + '_mergedcountstable_summary.txt', sep='\t')
    
    if generatePlots != 'off' and max(exptGroups.count().iloc[0]) > 1:
        printNow('-generating scatter plots of counts pre-merger')
    
        tempDataDict = {'library': libraryTable[sublibColumn],
                        'premerged counts': countsTable,
                       'counts': mergedCountsTable}

        for (phenotype, replicate), countsCols in exptGroups:
            if len(countsCols.columns) == 1:
                continue
            
            else:
                screen_analysis.premergedCountsScatterMatrix(tempDataDict, phenotype, replicate)

    if generatePlots != 'off':
        printNow('-generating sgRNA read count histograms')
    
        tempDataDict = {'library': libraryTable[sublibColumn],
                        'counts': mergedCountsTable}
                    
        for (phenotype, replicate), countsCol in mergedCountsTable.iteritems():
            screen_analysis.countsHistogram(tempDataDict, phenotype, replicate)
    
    #create pairs of columns for each comparison, filter to na, then generate sgRNA phenotype score
    printNow('Computing sgRNA phenotype scores')

    growthValueDict = {(tup[0],tup[1]):tup[2] for tup in exptParameters['growth_value_tuples']}
    phenotypeList = list(set(zip(*exptParameters['condition_tuples'])[0]))
    replicateList = sorted(list(set(zip(*exptParameters['counts_file_list'])[1])))

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
    
    if generatePlots  != 'off':
        tempDataDict = {'library': libraryTable[sublibColumn],
                        'counts': mergedCountsTable,
                        'phenotypes': pd.DataFrame(phenotypeScoreDict)}
                        
        printNow('-generating phenotype histograms and scatter plots')
        
        for (phenotype, condition1, condition2) in exptParameters['condition_tuples']:
            for replicate in replicateList:
                screen_analysis.countsScatter(tempDataDict, condition1, replicate, condition2, replicate, 
                    colorByPhenotype_condition = phenotype, colorByPhenotype_replicate = replicate)
                    
                screen_analysis.phenotypeHistogram(tempDataDict, phenotype, replicate)
                screen_analysis.sgRNAsPassingFilterHist(tempDataDict, phenotype, replicate)
    
    #scatterplot sgRNAs for all replicates, then average together and add columns to phenotype score table
    if len(replicateList) > 1:
        printNow('Averaging replicates')

        for phenotype in phenotypeList:
            repCols = pd.DataFrame({(phen,rep):col for (phen,rep), col in phenotypeScoreDict.iteritems() if phen == phenotype})
            phenotypeScoreDict[(phenotype,'ave_' + '_'.join(replicateList))] = repCols.mean(axis=1,skipna=False) #average nan and real to nan; otherwise this could lead to data points with just one rep informing results

    phenotypeTable = pd.DataFrame(phenotypeScoreDict)
    phenotypeTable.to_csv(outbase + '_phenotypetable.txt', sep='\t', tupleize_cols = False)

    if len(replicateList) > 1 and generatePlots != 'off':
        tempDataDict = {'library': libraryTable[sublibColumn],
                        'phenotypes': phenotypeTable}
                    
        printNow('-generating replicate phenotype histograms and scatter plots')
    
        for phenotype, phengroup in phenotypeTable.groupby(level=0, axis=1):
            for i, ((p, rep1), col1) in enumerate(phengroup.iteritems()):
                if rep1[:4] == 'ave_':
                    screen_analysis.phenotypeHistogram(tempDataDict, phenotype, rep1)
            
                for j, ((p, rep2), col2) in enumerate(phengroup.iteritems()):
                    if rep2[:4] == 'ave_' or j<=i:
                        continue
                    
                    else:
                        screen_analysis.phenotypeScatter(tempDataDict, phenotype, rep1, phenotype, rep2)                    
                

    #generate pseudogenes
    negTable = phenotypeTable.loc[libraryTable[sublibColumn].loc[:,'gene'] == 'negative_control',:]

    if exptParameters['generate_pseudogene_dist'] != 'off' and len(exptParameters['analyses']) > 0:
        print 'Generating a pseudogene distribution from negative controls'
        sys.stdout.flush()

        pseudoTableList = []
        pseudoLibTables = []
        negValues = negTable.values
        negColumns = negTable.columns

        if exptParameters['generate_pseudogene_dist'].lower() == 'manual':
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

        elif exptParameters['generate_pseudogene_dist'].lower() == 'auto':
            for pseudogene, (gene, group) in enumerate(libraryTable[sublibColumn].drop_duplicates(['gene','sequence']).groupby('gene')):
                if gene == 'negative_control':
                    continue 
                for transcript, (transcriptName, transcriptGroup) in enumerate(group.groupby('transcripts')):
                    randIndices = np.random.randint(0, len(negTable), len(transcriptGroup))
                    pseudoTable = negValues[randIndices,:]
                    pseudoIndex = ['pseudo_%d_%d_%d' % (pseudogene, transcript, i) for i in range(len(transcriptGroup))]
                    pseudoSeqs = ['seq_%d_%d_%d' % (pseudogene, transcript, i) for i in range(len(transcriptGroup))]
                    pseudoTableList.append(pd.DataFrame(pseudoTable,index=pseudoIndex,columns=negColumns))
                    pseudoLib = pd.DataFrame({'gene':['pseudo_%d'%pseudogene]*len(transcriptGroup),
                        'transcripts':['pseudo_transcript_%d'%transcript]*len(transcriptGroup),
                        'sequence':pseudoSeqs},index=pseudoIndex)
                    pseudoLibTables.append(pseudoLib)

        else:
            print 'generate_pseudogene_dist parameter not recognized, defaulting to off'

        phenotypeTable = phenotypeTable.append(pd.concat(pseudoTableList))
        libraryTableGeneAnalysis = libraryTable[sublibColumn].append(pd.concat(pseudoLibTables))
    else:
        libraryTableGeneAnalysis = libraryTable[sublibColumn]

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
    
    if generatePlots != 'off':
        if 'calculate_ave' in exptParameters['analyses'] and 'calculate_mw' in exptParameters['analyses']:
            tempDataDict = {'library': libraryTable[sublibColumn],
                            'gene scores': geneTableCollapsed if exptParameters['collapse_to_transcripts'] else geneTable}
                            
            for (phenotype, replicate), gtable in geneTableCollapsed.groupby(level=[0,1], axis=1):
                if len(replicateList) == 1 or replicate[:4] == 'ave_': #just plot averaged reps where available
                    screen_analysis.volcanoPlot(tempDataDict, phenotype, replicate, labelHits=True)

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


#compute phenotype scores for any given comparison of two conditions
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

    #compute neg control log2 enrichment
    if normToNegs == True:
        negCounts = combinedCountsPseudo.align(libraryTable[libraryTable['gene'] == 'negative_control'],axis=0,join='inner')[0]
        #print negCounts
    else:
        negCounts = combinedCountsPseudo
    neglog2e = negCounts.apply(calcLog2e, countsRatio=countsRatio, growthValue=1, wtLog2E=0, axis=1).median()
    #print neglog2e

    #compute phenotype scores
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

#apply gene scoring functions to pre-grouped tables of phenotypes
def applyGeneScoreFunction(groupedPhenotypeTable, negativeTable, analysis, analysisParamList):
    if analysis == 'calculate_ave':
        numToAverage = analysisParamList[0]
        if numToAverage <= 0:
            means = groupedPhenotypeTable.aggregate(np.mean)
            counts = groupedPhenotypeTable.count()
            result = pd.concat([means,counts],axis=1,keys=['average of all phenotypes','average of all phenotypes_sgRNAcount'])
        else:
            means = groupedPhenotypeTable.apply(lambda x: averageBestN(x, numToAverage))
            counts = groupedPhenotypeTable.count()
            result = pd.concat([means,counts],axis=1,keys=['average phenotype of strongest %d'%numToAverage, 'sgRNA count_avg'])
    elif analysis == 'calculate_mw':
        pvals = groupedPhenotypeTable.apply(lambda x: applyMW(x, negativeTable))
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

def averageBestN(group, numToAverage):
    return group.apply(lambda column: np.mean(sorted(column.dropna(),key=abs,reverse=True)[:numToAverage]) if len(column.dropna()) > 0 else np.nan)

def applyMW(group, negativeTable):
    if int(sp.__version__.split('.')[1]) >= 17: #implementation of the "alternative flag":
        return group.apply(lambda column: stats.mannwhitneyu(column.dropna().values, negativeTable[column.name].dropna().values, alternative = 'two-sided')[1] if len(column.dropna()) > 0 else np.nan)
    else:
        return group.apply(lambda column: stats.mannwhitneyu(column.dropna().values, negativeTable[column.name].dropna().values)[1] * 2 if len(column.dropna()) > 0 else np.nan) #pre v0.17 stats.mannwhitneyu is one-tailed!!


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



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate sgRNA- and gene-level phenotypes based on sequencing read counts, as specified by the experiment config file.')
    parser.add_argument('Config_File', help='Experiment config file specifying screen analysis settings (see accomapnying BLANK and DEMO files).')
    parser.add_argument('Library_File_Directory', help='Directory containing reference library tables and the library_config.txt file.')

    parser.add_argument('--plot_extension', default='png', help='Image extension for plot files, or \"off\". Default is png.')

    args = parser.parse_args()
    # print args

    processExperimentsFromConfig(args.Config_File, args.Library_File_Directory, args.plot_extension.lower())

