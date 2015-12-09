# Unintended_Consequences
Scripts and small data files for Unintended Consequences. In this project, we are comparing the patterns of genomic variation between standing variation in soybean cultivars, induced variation from fast neutron mutagenesis, and somaclonal variation from genetic transformation and tissue culture. Cultivar data comes from [Anderson et al. 2014](http://www.g3journal.org/content/4/7/1307.long), and fast neutron mutagenesis data comes from [Bolon et al. 2014](http://www.genetics.org/content/198/3/967). Transfomration data was collected for this study.

## Manuscript Status
*Submitted*

## Contents
- Data/
    - **FN.conf** - Configuration file for making a [Circos](http://circos.ca/) plot. Shows positions of induced base substitutions in FN-mutagenized lines.
    - **Mutation_Rates.txt**: Line-by-line generation, dosage, and number of substitutions 
    - **WPT.conf** - Configuration file for making a [Circos](http://circos.ca/) plot showing the positions of induced base substitutions in GTTC lines.
- Methods/
    - **Estimating_Mutation_Rate.pdf**: Supplemental file describing methods for estimating induced base substitution rate. Introduces terminology, assumptions, and equations used to calculate the rate of base substitution.
    - **Estimating_Mutation_Rate.tex**: TeX source file for the mutation rate estimation document.
    - **Spontaneous_Mutations.pdf**: PDF figure of the equations used in the supplemental methods document.
- Scripts/
    - **Mutation_Rate.R**: Calculates the rate of induced base substitution given the number of generations of inbreeding and the number of observed homozygous variants.
    - **ParseVCF.py**: Converts a VCF into a CSV that can be manipulated in R.
    - **Plot_Rates.R**: Makes a bar plot of rates of induced mutation by class.
    - **SNPHomozygousCalls.R** - Reads CSV from ``ParseVCF.py`` and calculates the number of homozygous mutations in each line.
    - **SNPcalling_Pipeline.sh** - Sequence analysis pipeline, will produce a VCF from a list of .FASTQ files.
