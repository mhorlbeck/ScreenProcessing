
# ScreenProcessing

##There are two steps to analyzing sequencing data from pooled screens using the ScreenProcessing pipeline
1. Convert raw sequencing files into library counts using bowtie alignment

    The script used for this step is **fastqgz_to_counts.py**, and requires a reference fasta file of the library used in the
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

(ScreenProcessing no longer uses Bowtie to align sequencing reads; if you want to use or fork from this functionality use an earlier version of the program)

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
```
usage: fastqgz_to_counts.py [-h]
                            [-u UNALIGNED_INDICES [UNALIGNED_INDICES ...]]
                            [-p PROCESSORS] [--trim_start TRIM_START]
                            [--trim_end TRIM_END] [--test]
                            Library_Fasta Out_File_Path Seq_File_Names
                            [Seq_File_Names ...]

Process raw sequencing data from screens to counts files in parallel

positional arguments:
  Library_Fasta         Fasta file of expected library reads.
  Out_File_Path         Directory where output files should be written.
  Seq_File_Names        Name(s) of sequencing file(s). Unix wildcards can be
                        used to select multiple files at once. The script will
                        search for all *.fastq.gz, *.fastq, and
                        *.fa(/fasta/fna) files with the given wildcard name.

optional arguments:
  -h, --help            show this help message and exit
  -p PROCESSORS, --processors PROCESSORS
  --trim_start TRIM_START
  --trim_end TRIM_END
  --test                Run the entire script on only the first 10000 reads of
                        each file. Be sure to delete or move all test files
                        before re-running script as they will not be
                        overwritten.
```

Let's try a test run--we'll write the command out but include the "--test" flag so we can catch any errors in the script that crop up along the way without waiting for the full dataset to be processed. Here's the command:

    python fastqgz_to_counts.py -p 6 --test --trim_start 1 --trim_end 35 bowtie_indices/CRISPRi_v1_human.trim_1_35.fa DEMO/output_folder_test sequencing_runs/140509_SR50/cholera/Sample_OC35_index*/*.fastq.gz
    
In order, the components of this command are:
* p 6 = we allocate up to 6 processors to running the script (this depends on your computer's capabilities and current load, and uses up to 1 processer per sequencing file so this won't )
* trim_start 1 = we trim off the first base of the sequencing read, because in this case first base of all of the reads were Gs and this "monotemplating" lowered the sequencing quality of the base
* trim_end 35 = we trim the read on the 3' end to the 35th base (as measured from the 5' end)

```
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index1/OC35_index1_ATCACG_L008_R1_001.fastq.gz
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index10/OC35_index10_TAGCTT_L008_R1_001.fastq.gz
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index12/OC35_index12_CTTGTA_L008_R1_001.fastq.gz
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index14/OC35_index14_AGTTCC_L008_R1_001.fastq.gz
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index3/OC35_index3_TTAGGC_L008_R1_001.fastq.gz
Processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index6/OC35_index6_GCCAAT_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index1/OC35_index1_ATCACG_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index6/OC35_index6_GCCAAT_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index12/OC35_index12_CTTGTA_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index10/OC35_index10_TAGCTT_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index3/OC35_index3_TTAGGC_L008_R1_001.fastq.gz
Done processing sequencing_runs/140509_SR50/cholera/Sample_OC35_index14/OC35_index14_AGTTCC_L008_R1_001.fastq.gz
DEMO/output_folder_test/count_files/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.26E+03 aligning (82.60%)
DEMO/output_folder_test/count_files/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.15E+03 aligning (81.47%)
DEMO/output_folder_test/count_files/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.13E+03 aligning (81.33%)
DEMO/output_folder_test/count_files/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.11E+03 aligning (81.14%)
DEMO/output_folder_test/count_files/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.11E+03 aligning (81.09%)
DEMO/output_folder_test/count_files/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	1.00E+04 reads	8.13E+03 aligning (81.30%)
```

The script ran successfully! The typical alignment rates we see for this step are in the 80-90% range, although if you just look at the first 10,000 reads (as we did using the --test flag) the percent aligning can be lower than this. 
Here's the output of the full run, as specified by just removing the --test flag:
    python fastqgz_to_counts.py -p 6 --trim_start 1 --trim_end 35 bowtie_indices/CRISPRi_v1_human.trim_1_35.fa DEMO/output_folder_full ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index*/*.fastq.gz
```
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index1/OC35_index1_ATCACG_L008_R1_001.fastq.gz
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index10/OC35_index10_TAGCTT_L008_R1_001.fastq.gz
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index12/OC35_index12_CTTGTA_L008_R1_001.fastq.gz
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index14/OC35_index14_AGTTCC_L008_R1_001.fastq.gz
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index3/OC35_index3_TTAGGC_L008_R1_001.fastq.gz
Processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index6/OC35_index6_GCCAAT_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index10/OC35_index10_TAGCTT_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index1/OC35_index1_ATCACG_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index6/OC35_index6_GCCAAT_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index14/OC35_index14_AGTTCC_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index12/OC35_index12_CTTGTA_L008_R1_001.fastq.gz
Done processing ../sequencing_runs/140509_SR50/cholera/Sample_OC35_index3/OC35_index3_TTAGGC_L008_R1_001.fastq.gz
DEMO/output_folder_full/count_files/OC35_index1_ATCACG_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.34E+07 reads	2.72E+07 aligning (81.31%)
DEMO/output_folder_full/count_files/OC35_index10_TAGCTT_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.35E+07 reads	2.71E+07 aligning (81.10%)
DEMO/output_folder_full/count_files/OC35_index12_CTTGTA_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.65E+07 reads	2.95E+07 aligning (80.85%)
DEMO/output_folder_full/count_files/OC35_index14_AGTTCC_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.46E+07 reads	2.80E+07 aligning (81.14%)
DEMO/output_folder_full/count_files/OC35_index3_TTAGGC_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.71E+07 reads	3.01E+07 aligning (81.04%)
DEMO/output_folder_full/count_files/OC35_index6_GCCAAT_L008_R1_001_CRISPRi_v1_human.trim_1_35.fa.counts:
	3.54E+07 reads	2.86E+07 aligning (80.77%)
Done processing all sequencing files
```

The unaligned_reads files represent the trimmed reads that did not match the library. These can be useful to diagnose any contaminants or systematic errors in your sequencing. These files can be quite large, and once the script is finished they aren't necessary, so once everything looks good you can delete those to save disk space. I've only left the counts_files in the DEMO/output_folder_full directory, since we'll need those for the next step. If you want to see what these files should look like, much smaller versions of them remain in DEMO/output_folder_test.

###Running process_experiments.py

This step turns the raw sgRNA counts into quantitative phenotypes for each sgRNA and gene in the experiment. Much of the analysis pipeline is drawn directly from the framework of Bassik and Kampmann et al Cell 2014 and Kampmann, Bassik, and Weissman PNAS 2014. 

To tell the script which sequencing files to use, you first need to fill out an experiment_config_file. A blank version is included with the scripts, along with a filled out example for the cholera toxin screens. Once this is filled out, run the script on the command line. The command line function is:
```
    python ScreenProcessing/process_experiments.py ScreenProcessing/experiment_config_file_DEMO.txt ScreenProcessing/library_tables/
    
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
```
The results are a set of tables, which can be opened in excel, or further analyzed using python or R. For interactive analysis, I strongly recommend the iPython Notebook (aka jupyter), which I also used to write this demo. The tables represent:
* the librarytable is the info for the sgRNAs
* the rawcountstable is simply the counts for each sgRNA from each sequencing file, while mergedcountsfile represents the counts after combining any experiments that may have been sequenced in multiple runs
* the phenotypetable is individual sgRNA phenotypes for replicates and averaged reps (gamma = growth, rho = toxin)
* the genetable is gene scores/p-values for each TSS
* the genetable_collapsed is gene scores for each gene (using the TSS with the best p-value)

If you followed along using just the 6 sequencing files I specified at the start, your results won't exactly match the published results using the full complement of sequencing reads, but should match the tables included here.

Functions and examples for downstream analysis are in progress.
