from ConfigParser import SafeConfigParser
import os
import fnmatch

#Parse and validate the input of an experiment config file
#output a dict with all of the parameters needed to process experiments
def parseExptConfig(configFile, librariesToSublibrariesDict):
    parser = SafeConfigParser()
    results = parser.read(configFile)
    if len(results) == 0:
        return None, 1, 'Experiment config file not found'
    
    #output variables
    paramDict = dict()
    exitStatus = 0
    warningString = ''

    ##check all sections
    expectedSections = set(['experiment_settings',
         'library_settings',
         'counts_files',
         'filter_settings',
         'sgrna_analysis',
         'growth_values',
         'gene_analysis'])
         
    parsedSections = set(parser.sections())
    
    if len(expectedSections) != len(parsedSections) and len(expectedSections) != len(expectedSections.intersection(parsedSections)):
        return paramDict, 1, 'Config file does not have all required sections or has extraneous sections!\nExpected:' + ','.join(expectedSections) + '\nFound:' + ','.join(parsedSections)

    ##experiment settings
    if parser.has_option('experiment_settings','output_folder'):
        paramDict['output_folder'] = parser.get('experiment_settings','output_folder') #ways to check this is a valid path?
    else:
        warningString += 'No output folder specified, defaulting to current directory\n.'
        paramDict['output_folder'] = os.curdir()
        
    if parser.has_option('experiment_settings','experiment_name'):
        paramDict['experiment_name'] = parser.get('experiment_settings','experiment_name')
    else:
        warningString += 'No experiment name specified, defaulting to \'placeholder_expt_name\'\n.'
        paramDict['experiment_name'] = 'placeholder_expt_name'
        
        
    ##library settings
    libraryDict = librariesToSublibrariesDict
    if parser.has_option('library_settings','library'):
        parsedLibrary = parser.get('library_settings','library')
        
        if parsedLibrary.lower() in libraryDict:
            paramDict['library'] = parsedLibrary.lower()
        else:
            warningString += 'Library name \"%s\" not recognized\n' % parsedLibrary
            exitStatus += 1
            
    else:
        warningString += 'No library specified\n'
        exitStatus += 1
        parsedLibrary = ''
      
    if 'library' in paramDict:
        if parser.has_option('library_settings','sublibraries'):
            parsedSubList = parser.get('library_settings','sublibraries').strip().split('\n')
            
            paramDict['sublibraries'] = []
            
            for sub in parsedSubList:
                sub = sub.lower()
                if sub in libraryDict[paramDict['library']]:
                    paramDict['sublibraries'].append(sub)
                    
                else:
                    warningString += 'Sublibrary %s not recognized\n' % sub
                    
        else:
            paramDict['sublibraries'] = libraryDict[paramDict['library']]
            
    
    ##counts files
    if parser.has_option('counts_files','counts_file_string'):
        countsFileString = parser.get('counts_files','counts_file_string').strip()
        
        paramDict['counts_file_list'] = []
        
        for stringLine in countsFileString.split('\n'):
            stringLine = stringLine.strip()
            
            if len(stringLine.split(':')) != 2 or len(stringLine.split('|')) != 2:
                warningString += 'counts file entry could not be parsed: ' + stringLine + '\n'
                exitStatus += 1
            
            else:
                parsedPath = stringLine.split(':')[0]
                
                if os.path.isfile(parsedPath) == False:
                    warningString += 'Counts file not found: ' + parsedPath + '\n'
                    exitStatus += 1
                    
                condition, replicate = stringLine.split(':')[1].split('|')
                
                paramDict['counts_file_list'].append((condition, replicate,parsedPath))
            
    else:
        warningString += 'No counts files entered\n'
        exitStatus += 1
        
    
    ##filter settings
    filterOptions = ['either','both']
    if parser.has_option('filter_settings','filter_type') and parser.get('filter_settings','filter_type').lower() in filterOptions:
        paramDict['filter_type'] = parser.get('filter_settings','filter_type').lower()
    else:
        warningString += 'Filter type not set or not recognized, defaulting to \'either\'\n'
        paramDict['filter_type'] = 'either'
        
    if parser.has_option('filter_settings','minimum_reads'):
        try:
            paramDict['minimum_reads'] = parser.getint('filter_settings','minimum_reads')
        except ValueError:
            warningString += 'Minimum read value not an integer, defaulting to 0\n'  #recommended value is 50 but seems arbitrary to default to that
            paramDict['minimum_reads'] = 0
    else:
        warningString += 'Minimum read value not found, defaulting to 0\n'  #recommended value is 50 but seems arbitrary to default to that
        paramDict['minimum_reads'] = 0


    ##sgRNA Analysis
    if parser.has_option('sgrna_analysis','condition_string'):
        conditionString = parser.get('sgrna_analysis','condition_string').strip()

        paramDict['condition_tuples'] = []
        
        if 'counts_file_list' in paramDict:
            expectedConditions = set(zip(*paramDict['counts_file_list'])[0])
        else:
            expectedConditions = []

        enteredConditions = set()
        
        for conditionStringLine in conditionString.split('\n'):
            conditionStringLine = conditionStringLine.strip()

            if len(conditionStringLine.split(':')) != 3:
                warningString += 'Phenotype condition line not understood: ' + conditionStringLine + '\n'
                exitStatus += 1
            else:
                phenotype, condition1, condition2 = conditionStringLine.split(':')

                if condition1 not in expectedConditions or condition2 not in expectedConditions:
                    warningString += 'One of the conditions entered does not correspond to a counts file: ' + conditionStringLine + '\n'
                    exitStatus += 1
                else:
                    paramDict['condition_tuples'].append((phenotype,condition1,condition2))
                    enteredConditions.add(condition1)
                    enteredConditions.add(condition2)

        if len(paramDict['condition_tuples']) == 0:
            warningString += 'No phenotype score/condition pairs found\n'
            exitStatus += 1
            
        unusedConditions = list(expectedConditions - enteredConditions)
        if len(unusedConditions) > 0:
            warningString += 'Some conditions assigned to counts files will not be incorporated in sgRNA analysis:\n' \
                + ','.join(unusedConditions) + '\n'
        
    else:
        warningString += 'No phenotype score/condition pairs entered\n'
        exitStatus += 1
    
    
    pseudocountOptions = ['zeros only','all values','filter out']
    if parser.has_option('sgrna_analysis','pseudocount_behavior') and parser.get('sgrna_analysis','pseudocount_behavior').lower() in pseudocountOptions:
        paramDict['pseudocount_behavior'] = parser.get('sgrna_analysis','pseudocount_behavior').lower()
    else:
        warningString += 'Pseudocount behavior not set or not recognized, defaulting to \'zeros only\'\n'
        paramDict['pseudocount_behavior'] = 'zeros only'

    if parser.has_option('sgrna_analysis','pseudocount'):
        try:
            paramDict['pseudocount'] = parser.getfloat('sgrna_analysis','pseudocount')
        except ValueError:
            warningString += 'Pseudocount value not an number, defaulting to 0.1\n'
            paramDict['pseudocount'] = 0.1
    else:
        warningString += 'Pseudocount value not found, defaulting to 0.1\n'
        paramDict['pseudocount'] = 0.1

    
    ##Growth Values
    if parser.has_option('growth_values','growth_value_string') and len(parser.get('growth_values','growth_value_string').strip()) != 0:
        growthValueString = parser.get('growth_values','growth_value_string').strip()

        if 'condition_tuples' in paramDict and 'counts_file_list' in paramDict:
            expectedComparisons = set(zip(*paramDict['condition_tuples'])[0])
            expectedReplicates = set(zip(*paramDict['counts_file_list'])[1])

            expectedTupleList = []

            for comp in expectedComparisons:
                for rep in expectedReplicates:
                    expectedTupleList.append((comp,rep))
        else:
            expectedTupleList = []
             
        enteredTupleList = []
        growthValueTuples = []
        
        for growthValueLine in growthValueString.split('\n'):
            growthValueLine = growthValueLine.strip()
            
            linesplit = growthValueLine.split(':')
            
            if len(linesplit) != 3:
                warningString += 'Growth value line not understood: ' + growthValueLine + '\n'
                exitStatus += 1
                continue
                
            comparison = linesplit[0]
            replicate = linesplit[1]
            
            try:
                growthVal = float(linesplit[2])
            except ValueError:
                warningString += 'Growth value not a number: ' + growthValueLine + '\n'
                exitStatus += 1
                continue

            curTup = (comparison,replicate)
            if curTup in expectedTupleList:
                if curTup not in enteredTupleList:
                    enteredTupleList.append(curTup)
                    growthValueTuples.append((comparison,replicate,growthVal))
                    
                else:
                    warningString += ':'.join(curTup) + ' has multiple growth values entered\n'
                    exitStatus += 1
            else:
                warningString += ':'.join(curTup) + ' was not expected given the specified counts file assignments and sgRNA phenotypes\n'
                exitStatus += 1
        
        #because we enforced no duplicates or unexpected values these should match up unless there were values not entered
        #require all growth values to be explictly entered if some were
        if len(enteredTupleList) != len(expectedTupleList):
            warningString += 'Growth values were not entered for all expected comparisons/replicates. Expected: ' + \
                ','.join([':'.join(tup) for tup in expectedTupleList]) + '\nEntered: ' + \
                ','.join([':'.join(tup) for tup in enteredTupleList]) + '\n'
            exitStatus += 1
        else:
            paramDict['growth_value_tuples'] = growthValueTuples
            
    else:
        warningString += 'No growth values--all phenotypes will be reported as log2enrichments\n'
        
        paramDict['growth_value_tuples'] = []
        
        if 'condition_tuples' in paramDict and 'counts_file_list' in paramDict:
            expectedComparisons = set(zip(*paramDict['condition_tuples'])[0])
            expectedReplicates = set(zip(*paramDict['counts_file_list'])[1])

            for comp in expectedComparisons:
                for rep in expectedReplicates:
                    paramDict['growth_value_tuples'].append((comp,rep,1))
    
    ##Gene Analysis
    if parser.has_option('gene_analysis','collapse_to_transcripts'):
        try:
            paramDict['collapse_to_transcripts'] = parser.getboolean('gene_analysis','collapse_to_transcripts')
        except ValueError:
            warningString += 'Collapse to transcripts entry not a recognized boolean value\n'
            exitStatus += 1
    else:
        paramDict['collapse_to_transcripts'] = True
        warningString += 'Collapse to transcripts defaulting to True\n'


    #pseudogene parameters
    if parser.has_option('gene_analysis','generate_pseudogene_dist'):
        paramDict['generate_pseudogene_dist'] = parser.get('gene_analysis','generate_pseudogene_dist').lower()

        if paramDict['generate_pseudogene_dist'] not in ['auto','manual','off']:
            warningString += 'Generate pseudogene dist entry not a recognized option\n'
            exitStatus += 1
    else:
        paramDict['generate_pseudogene_dist'] = False
        warningString += 'Generate pseudogene dist defaulting to False\n'

    if 'generate_pseudogene_dist' in paramDict and paramDict['generate_pseudogene_dist'] == 'manual':
        if parser.has_option('gene_analysis','pseudogene_size'):
            try:
                paramDict['pseudogene_size'] = parser.getint('gene_analysis','pseudogene_size')
            except ValueError:
                warningString += 'Pseudogene size entry not a recognized integer value\n'
                exitStatus += 1
        else:
            warningString += 'No pseudogene size provided\n'
            exitStatus += 1

        if parser.has_option('gene_analysis','num_pseudogenes'):
            try:
                paramDict['num_pseudogenes'] = parser.getint('gene_analysis','num_pseudogenes')
            except ValueError:
                warningString += 'Pseudogene number entry not a recognized integer value\n'
                exitStatus += 1
        else:
            warningString += 'No pseudogene size provided\n'

    #list possible analyses in param dict as dictionary with keys = analysis and values = analysis-specific params
    
    paramDict['analyses'] = dict()

    #analyze by average of best n
    if parser.has_option('gene_analysis','calculate_ave'):
        try:
            if parser.getboolean('gene_analysis','calculate_ave') == True:
                paramDict['analyses']['calculate_ave'] = []
        except ValueError:
            warningString += 'Calculate ave entry not a recognized boolean value\n'
            exitStatus += 1

        if 'calculate_ave' in paramDict['analyses']:
            if parser.has_option('gene_analysis','best_n'):
                try:
                    paramDict['analyses']['calculate_ave'].append(parser.getint('gene_analysis','best_n'))
                except ValueError:
                    warningString += 'Best_n entry not a recognized integer value\n'
                    exitStatus += 1
            else:
                warningString += 'No best_n value provided for average analysis function\n'
                exitStatus += 1
    else:
        warningString += 'Best n average analysis not specified, defaulting to False\n'
    
    #analyze by Mann-Whitney
    if parser.has_option('gene_analysis','calculate_mw'):
        try:
            if parser.getboolean('gene_analysis','calculate_mw') == True:
                paramDict['analyses']['calculate_mw'] = []
        except ValueError:
            warningString += 'Calculate Mann-Whitney entry not a recognized boolean value\n'
            exitStatus += 1

    #analyze by K-S, skipping for now

    #analyze by nth best sgRNA
    if parser.has_option('gene_analysis','calculate_nth'):
        try:
            if parser.getboolean('gene_analysis','calculate_nth') == True:
                paramDict['analyses']['calculate_nth'] = []
        except ValueError:
            warningString += 'Calculate best Nth sgRNA entry not a recognized boolean value\n'
            exitStatus += 1

        if 'calculate_nth' in paramDict['analyses']:
            if parser.has_option('gene_analysis','nth'):
                try:
                    paramDict['analyses']['calculate_nth'].append(parser.getint('gene_analysis','nth'))
                except ValueError:
                    warningString += 'Nth best sgRNA entry not a recognized integer value\n'
                    exitStatus += 1
            else:
                warningString += 'No Nth best value provided for that analysis function\n'
                exitStatus += 1
    else:
        warningString += 'Nth best sgRNA analysis not specified, defaulting to False\n'


    if len(paramDict['analyses']) == 0:
        warningString += 'No analyses selected to compute gene scores\n' #should this raise exitStatus?

    return paramDict, exitStatus, warningString
    
#Parse the library configuration file to get the available libraries, sublibraries, and corresponding library table files
def parseLibraryConfig(libConfigFile):
    parser = SafeConfigParser()
    result = parser.read(libConfigFile)
    if len(result) == 0:
        raise ValueError('Library config file not found')
    
    librariesToSublibraries = dict()
    librariesToTables = dict()
    for library in parser.sections():
        tableFile = parser.get(library, 'filename').strip()
        librariesToTables[library.lower()] = tableFile

        sublibraryList = parser.get(library, 'sublibraries').strip().split('\n')
        librariesToSublibraries[library.lower()] = [sub.strip().lower() for sub in sublibraryList]

    if len(librariesToTables) == 0:
        raise ValueError('Library config file empty')
    
    return librariesToSublibraries, librariesToTables
