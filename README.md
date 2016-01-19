
# ScreenProcessing

##There are two steps to analyzing sequencing data from pooled screens using the ScreenProcessing pipeline
1. Convert raw sequencing files into library counts using bowtie alignment

    The script used for this step is **fastqgz_to_counts.py**, and requires a bowtie index of the library used in the
    experiment. Indices for several published libraries are included here, and more can be generated upon request.
    <p>
2. Generate sgRNA phenotype scores, gene-level scores, and gene-level p-values

    This step relies on the script **process_experiments.py**. This script requires you to first fill out a configuration file which allows you to:
    * Assign the read counts generated for each sequencing file to the appropriate sample condition
    * Specify which conditions you want to compare in order to generate phenotypes
    * Set several data quality filtering and normalization parameters
    * Specify how to score genes

###Dependencies
* Python v2.7
* Biopython
* Scipy/Numpy/Pandas/Matplotlib

* Bowtie1 in Unix PATH

#ScreenProcessing Demo

*If you want to follow along at home, you can access the data here:* [data link](https://ucsf.box.com/s/gq1lsrrl1eaz9ur0j5zc6ww2ebfx24zn)

The files correspond to the full complement of sequencing data used for the cell growth and cholera toxin sensitivity CRISPRi screens published in Gilbert and Horlbeck et al., Cell 2014 (we sequenced deeply to test out the system!) but for the purposes of the demo you can just access
* OC35_index1_ATCACG_L008_R1_001.fastq.gz
* OC35_index10_TAGCTT_L008_R1_001.fastq.gz
* OC35_index12_CTTGTA_L008_R1_001.fastq.gz
* OC35_index14_AGTTCC_L008_R1_001.fastq.gz
* OC35_index3_TTAGGC_L008_R1_001.fastq.gz
* OC35_index6_GCCAAT_L008_R1_001.fastq.gz

###Running fastqgz_to_counts.py

Open the terminal and go to the folder containing the script. You should be able to type:

    python fastqgz_to_counts.py -h

and see the usage information for running the script:

    usage: fastqgz_to_counts.py [-h]
                                [-u UNALIGNED_INDICES [UNALIGNED_INDICES ...]]
                                [-p PROCESSORS] [--trim_start TRIM_START]
                                [--trim_end TRIM_END] [--test]
                                Bowtie_Index Out_File_Path Seq_File_Names
                                [Seq_File_Names ...]

    Process raw sequencing data from screens to counts files in parallel

    positional arguments:
      Bowtie_Index          Bowtie index root name (without .1.ebwt).
      Out_File_Path
      Seq_File_Names        Name(s) of sequencing file(s). Unix wildcards can be
                            used to select multiple files at once. The script will
                            search for all *.fastq.gz, *.fastq, and
                            *.fa(/fasta/fna) files with the given wildcard name.

    optional arguments:
      -h, --help            show this help message and exit
      -u UNALIGNED_INDICES [UNALIGNED_INDICES ...], --Unaligned_Indices UNALIGNED_INDICES [UNALIGNED_INDICES ...]
                            Bowtie indices to test unaligned reads for possible
                            cross-contaminations.
      -p PROCESSORS, --processors PROCESSORS
      --trim_start TRIM_START
      --trim_end TRIM_END
      --test                Run the entire script on only the first 100 reads of
                            each file. Be sure to delete or move all test files
                            before re-running script as they will not be
                            overwritten.

    

Let's try a test run--we'll write the command out but include the "--test" flag so we can catch any errors in the script that crop up along the way without waiting for the full dataset to be processed. Here's the command:

    python fastqgz_to_counts.py --test -p 12 --trim_start 1 --trim_end 35 bowtie_indices/CRISPRi_v1_human DEMO/output_folder_test sequencing_runs/140509_SR50/cholera/Sample_OC35_index*/*.fastq.gz
    
In order, the components of this command are:
* p 12 = we allocate up to 12 processors to running the script (this depends on your computer's capabilities and current load)
* trim_start 1 = we trim off the first base of the sequencing read, because in this case first base of all of the reads were Gs and this "monotemplating" lowered the sequencing quality of the base
* trim_end 35 = we trim the read on the 3' end to the 35th base (as measured from the 5' end)

    
    Parsing and trimming sequence files...
    output_folder/trimmed_fastas/OC35_index10_TAGCTT_L008_R1_001_trim.fa:	100 reads
    output_folder/trimmed_fastas/OC35_index12_CTTGTA_L008_R1_001_trim.fa:	100 reads
    output_folder/trimmed_fastas/OC35_index14_AGTTCC_L008_R1_001_trim.fa:	100 reads
    output_folder/trimmed_fastas/OC35_index1_ATCACG_L008_R1_001_trim.fa:	100 reads
    output_folder/trimmed_fastas/OC35_index3_TTAGGC_L008_R1_001_trim.fa:	100 reads
    output_folder/trimmed_fastas/OC35_index6_GCCAAT_L008_R1_001_trim.fa:	100 reads
    Done parsing and trimming sequence files
    Mapping all sequencing runs to index: bowtie_indices/CRISPRi_v1_human
    Starting output_folder/bowtiemaps/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 83 (83.00%)
    # reads that failed to align: 17 (17.00%)
    Reported 86 alignments to 1 output stream(s)
    Starting output_folder/bowtiemaps/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 78 (78.00%)
    # reads that failed to align: 22 (22.00%)
    Reported 80 alignments to 1 output stream(s)
    Starting output_folder/bowtiemaps/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 75 (75.00%)
    # reads that failed to align: 25 (25.00%)
    Reported 76 alignments to 1 output stream(s)
    Starting output_folder/bowtiemaps/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 84 (84.00%)
    # reads that failed to align: 16 (16.00%)
    Reported 84 alignments to 1 output stream(s)
    Starting output_folder/bowtiemaps/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 73 (73.00%)
    # reads that failed to align: 27 (27.00%)
    Reported 75 alignments to 1 output stream(s)
    Starting output_folder/bowtiemaps/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 100
    # reads with at least one reported alignment: 74 (74.00%)
    # reads that failed to align: 26 (26.00%)
    Reported 77 alignments to 1 output stream(s)
    Done bowtie mapping
    Counting mapped reads...
    output_folder/count_files/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.counts: 	8.600000e+01 total mappings	2.000000e+00 ambiguous reads
    output_folder/count_files/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.counts: 	8.000000e+01 total mappings	2.000000e+00 ambiguous reads
    output_folder/count_files/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.counts: 	7.600000e+01 total mappings	1.000000e+00 ambiguous reads
    output_folder/count_files/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.counts: 	8.400000e+01 total mappings	0.000000e+00 ambiguous reads
    output_folder/count_files/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.counts: 	7.500000e+01 total mappings	2.000000e+00 ambiguous reads
    output_folder/count_files/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.counts: 	7.700000e+01 total mappings	3.000000e+00 ambiguous reads
    Done counting mapped reads

The script ran successfully! The typical alignment rates we see for this step are in the 80-90% range, although if you just look at the first 100 reads (as we did using the --test flag) the percent aligning is often much lower. 
Here's the output of the full run, as specified by just removing the --test flag:
    python fastqgz_to_counts.py -p 12 --trim_start 1 --trim_end 35 bowtie_indices/CRISPRi_v1_human DEMO/output_folder_full sequencing_runs/140509_SR50/cholera/Sample_OC35_index*/*.fastq.gz

    Parsing and trimming sequence files...
    DEMO/output_folder_full/trimmed_fastas/OC35_index10_TAGCTT_L008_R1_001_trim.fa: 33450053 reads
    DEMO/output_folder_full/trimmed_fastas/OC35_index12_CTTGTA_L008_R1_001_trim.fa: 36513934 reads
    DEMO/output_folder_full/trimmed_fastas/OC35_index14_AGTTCC_L008_R1_001_trim.fa: 34556020 reads
    DEMO/output_folder_full/trimmed_fastas/OC35_index1_ATCACG_L008_R1_001_trim.fa:  33420630 reads
    DEMO/output_folder_full/trimmed_fastas/OC35_index3_TTAGGC_L008_R1_001_trim.fa:  37135629 reads
    DEMO/output_folder_full/trimmed_fastas/OC35_index6_GCCAAT_L008_R1_001_trim.fa:  35431784 reads
    Done parsing and trimming sequence files
    Mapping all sequencing runs to index: bowtie_indices/CRISPRi_v1_human
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 33450053
    # reads with at least one reported alignment: 27318601 (81.67%)
    # reads that failed to align: 6131452 (18.33%)
    Reported 45017835 alignments to 1 output stream(s)
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 36513934
    # reads with at least one reported alignment: 29729223 (81.42%)
    # reads that failed to align: 6784711 (18.58%)
    Reported 51860863 alignments to 1 output stream(s)
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 34556020
    # reads with at least one reported alignment: 28231126 (81.70%)
    # reads that failed to align: 6324894 (18.30%)
    Reported 46076448 alignments to 1 output stream(s)
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 33420630
    # reads with at least one reported alignment: 27360011 (81.87%)
    # reads that failed to align: 6060619 (18.13%)
    Reported 45735197 alignments to 1 output stream(s)
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 37135629
    # reads with at least one reported alignment: 30304183 (81.60%)
    # reads that failed to align: 6831446 (18.40%)
    Reported 52683811 alignments to 1 output stream(s)
    Starting DEMO/output_folder_full/bowtiemaps/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.map ...
    # reads processed: 35431784
    # reads with at least one reported alignment: 28819719 (81.34%)
    # reads that failed to align: 6612065 (18.66%)
    Reported 49618732 alignments to 1 output stream(s)
    Done bowtie mapping
    Counting mapped reads...
    DEMO/output_folder_full/count_files/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.counts:    4.501784e+07 total mappings     1.160882e+06 ambiguous reads
    DEMO/output_folder_full/count_files/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.counts:    5.186086e+07 total mappings     1.274742e+06 ambiguous reads
    DEMO/output_folder_full/count_files/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.counts:    4.607645e+07 total mappings     1.030809e+06 ambiguous reads
    DEMO/output_folder_full/count_files/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.counts:     4.573520e+07 total mappings     9.964510e+05 ambiguous reads
    DEMO/output_folder_full/count_files/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.counts:     5.268381e+07 total mappings     1.136886e+06 ambiguous reads
    DEMO/output_folder_full/count_files/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.counts:     4.961873e+07 total mappings     1.098383e+06 ambiguous reads
    Done counting mapped reads


The trimmed_fastas and bowtie_maps files are quite large, and once the script is finished they aren't necessary, so once everything looks good feel free to delete those to save disk space. I've only left the counts_files in the DEMO/output_folder_full directory, since we'll need those for the next step. If you want to see what the intermediate files should look like, much smaller versions of them remain in DEMO/output_folder_test.

###Running process_experiments.py

This step turns the raw sgRNA counts into quantitative phenotypes for each sgRNA and gene in the experiment. Much of the analysis pipeline is drawn directly from the framework of Bassik and Kampmann et al Cell 2014 and Kampmann, Bassik, and Weissman PNAS 2014. 

To tell the script which sequencing files to use, you first need to fill out an experiment_config_file. A blank version is included with the scripts, along with a filled out example for the cholera toxin screens. Once this is filled out, run the script using an interactive python interpreter (a command line function is in progress). I strongly recommend the iPython Notebook (aka jupyter), which I also used to write this demo.


    %run ScreenProcessing/process_experiments.py #load the script into the interpreter
    #run the command to process experiments, giving the path the the config file you wrote and the directory where the reference tables are kept
    processExperimentsFromConfig('ScreenProcessing/experiment_config_file_DEMO.txt','ScreenProcessing/library_tables/')

    
    Accessing library information
    Loading counts data
    Merging experiment counts split across lanes/indexes
    Generating scatter plots of counts pre-merger and histograms of counts post-merger
    Computing sgRNA phenotype scores
    Plotting and averaging replicates
    Generating a pseudogene distribution from negative controls
    Computing gene scores
    --calculate_ave
    --calculate_mw
    Collapsing transcript scores to gene scores
    Done!


    /usr/local/lib/python2.7/dist-packages/numpy/core/_methods.py:59: RuntimeWarning: Mean of empty slice.
      warnings.warn("Mean of empty slice.", RuntimeWarning)


The results are a set of tables, which can be opened in excel, or further analyzed using python or R. The tables represent:
* the librarytable is the info for the sgRNAs
* the rawcountstable is simply the counts for each sgRNA from each sequencing file, while mergedcountsfile represents the counts after combining any experiments that may have been sequenced in multiple runs
* the phenotypetable is individual sgRNA phenotypes for replicates and averaged reps (gamma = growth, rho = toxin)
* the genetable is gene scores/p-values for each TSS
* the genetable_collapsed is gene scores for each gene (using the TSS with the best p-value)

If you followed along using just the 6 sequencing files I specified at the start, your results won't exactly match the published results using the full complement of sequencing reads, but should match the tables included here.

Functions and examples for downstream analysis are in progress.
