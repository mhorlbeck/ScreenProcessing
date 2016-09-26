
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

    This script also generates a set of standard graphs using **screen_analysis.py**

3. [Optional] Generate graphs interactively using **screen_analysis.py**

###Dependencies
* Python v2.7
* Biopython
* Scipy/Numpy/Pandas/Matplotlib
* iPython or iPython Notebook recommended for interactive graph plotting

(ScreenProcessing no longer uses Bowtie to align sequencing reads; if you want to use or fork from this functionality use an earlier version of the program)

#ScreenProcessing Demo
A PDF slideshow with a step-by-step tutorial of screen analysis using the data files included in the Demo folder can found here: [ScreenProcessing Demo](ScreenProcessing_tutorial.pdf)

*The demo files represent a tiny slice of the full sequencing dataset to speed up the download and demo scripts. The full complement of sequencing data used for the cell growth and cholera toxin sensitivity CRISPRi screens published in Gilbert and Horlbeck et al., Cell 2014 can be accessed here:* [data link](https://ucsf.box.com/s/gq1lsrrl1eaz9ur0j5zc6ww2ebfx24zn)
