#######First thing that you need to do it to combine your forward read runs and reverse read runs using
#######the UNIX cat function. Make sure you preserve the adapter sequence in your final file name
####This is an example pipe line for one data set to generate SNPs in soybean
#go the diretory of where the Data is stored

#   Written by Jean-Michel Michno
cd /home/stuparr/shared/Projects/NGS/UMGC_Archive/141126_D00635_0029_AC57YDANXX/Project_Stupar_Project_022

#use the can function using a space inbetween each file you want to combine
#> write your cat  output to a .fastq file in your specified directory
cat \
WPT391-1-6_CGTACG_L001_R1_001.fastq \
WPT391-1-6_CGTACG_L002_R1_001.fastq \
WPT391-1-6_CGTACG_L003_R1_001.fastq \
WPT391-1-6_CGTACG_L004_R1_001.fastq \
WPT391-1-6_CGTACG_L005_R1_001.fastq \
WPT391-1-6_CGTACG_L006_R1_001.fastq \
WPT391-1-6_CGTACG_L007_R1_001.fastq \
WPT391-1-6_CGTACG_L008_R1_001.fastq \
> /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1.fastq

# do the same for your reverse reads
cat \
WPT391-1-6_CGTACG_L001_R2_001.fastq \
WPT391-1-6_CGTACG_L002_R2_001.fastq \
WPT391-1-6_CGTACG_L003_R2_001.fastq \
WPT391-1-6_CGTACG_L004_R2_001.fastq \
WPT391-1-6_CGTACG_L005_R2_001.fastq \
WPT391-1-6_CGTACG_L006_R2_001.fastq \
WPT391-1-6_CGTACG_L007_R2_001.fastq \
WPT391-1-6_CGTACG_L008_R2_001.fastq \
> /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2.fastq

#########################Clean your Reads
#go to my scratch directory
 cd /home/stuparr/michnoj0/Scratch/WPT391-1-6


#load and run fastqc module to check the quality of the reads
module load fastqc
fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1.fastq

#load and run cutadapt
#-b is your uniquie adapter with the 6bp line being the 6bp listed in your fastq file
# the second -b is the universal adater
#-f fastq is to specify a fastq input
#-m is your minimum read length allowed
#-quality-base is your phred score
#-o is you output file <space> then your input file
#
module load cutadapt
cutadapt \
-b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\
CGTACG\
ATCTCGTATGCCGTCTTCTGCTTG \
-b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
-f fastq \
-m 30 \
--quality-base=33 \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt.fastq \
/home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1.fastq \
> ~/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt.log
#run fastq to see how the quality has changed
fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt.fastq


#load and run fastx for low complexity sequences
module load fastx_toolkit
#-Q is for phred score
#-i is input
#-o is output
fastx_artifacts_filter -Q33 -v \
-i /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt.fastq \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity.fastq

#run fastq again to make sure you are trimming reads
fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity.fastq


#-Q is phred score
#-v reports trimming info
#-t is your minumum accepted phred score
#-l is your minimum accepted length
fastq_quality_trimmer -Q  33 -v -t 20 -l 30 \
-i /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity.fastq \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity_qual_trim.fastq 

#run a check again
fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity_qual_trim.fastq 


########### run the same thing but for the reverse reads

#go to my scratch directory
 cd /home/stuparr/michnoj0/Scratch/WPT391-1-6


#load and run fastqc module
module load fastqc
fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2.fastq


#load and run cutadapt
module load cutadapt
cutadapt \
-b GATCGGAAGAGCACACGTCTGAACTCCAGTCAC\
CGTACG\
ATCTCGTATGCCGTCTTCTGCTTG \
-b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
-f fastq \
-m 30 \
--quality-base=33 \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt.fastq \
/home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2.fastq \
> ~/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt.log

fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt.fastq


#load and run fastx for low complexity sequences
module load fastx_toolkit
fastx_artifacts_filter -Q33 -v \
-i /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt.fastq \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity.fastq

fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity.fastq



fastq_quality_trimmer -Q  33 -v -t 20 -l 30 \
-i /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity.fastq \
-o /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity_qual_trim.fastq 


fastqc -f fastq /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity_qual_trim.fastq 


##################################Pair and Sync your reads
###Sync your reads from the F and R so they can be properly alligned
module load riss_util

cd /home/stuparr/michnoj0/Scratch/WPT391-1-6

resync.pl \
WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity_qual_trim.fastq \
WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity_qual_trim.fastq

#your output will be <input>.out


#this tool will check to see that your reads are synced properly
module load pe-sync

pe-sync-2-files.pl \
WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity_qual_trim.fastq.out \
WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity_qual_trim.fastq.out 

################################ map your reads
#go to the directory where the fastq files will be located
 cd ~/Scratch/WPT391-1-6/

 #load modules
module load bowtie/2.2.4
module load samtools/0.1.18

#run Bowtie to do the end to end alignment
#-q is for fastq file
#-x the base name of the index of the reference genome
#-rg-id sets the read group ID
#--phred33 is to set the phred score of 33
#-p is to set how may processors/threads do you want to use
# -1 and -2 are the forward and reverse reads fastq files
#-S is what you want the output file to be called

bowtie2 \
-q \
-x /home/stuparr/shared/References/Gmax.a1.v1.1/assembly/Gmax_189 \
--rg-id WPT391-1-6 \
--phred33 \
-p 8 \
-1 /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR1_cutadapt_complexity_qual_trim.fastq.out \
-2 /home/stuparr/michnoj0/Scratch/WPT391-1-6/WPT391-1-6_CGTACG_MergedR2_cutadapt_complexity_qual_trim.fastq.out \
-S WPT391-1-6.sam

#convert the sam file to a bam file
samtools view -bS WPT391-1-6.sam > WPT391-1-6.bam

#sort the bam file
samtools sort WPT391-1-6.bam WPT391-1-6.sorted

#index the bam so you can view in IGV
samtools index WPT391-1-6.sorted

#You can now use the sorted bam and the .bai file together for view in IGV
















#then use mpileup to call variants.....

########################################
cd /home/stuparr/michnoj0/Scratch/WPT391-1-6/

 #load modules
module load samtools/0.1.18
module load bedtools/2.17.0

#use bedtools to get the coverage of the genome
#-dz Report the depth at each genome position with 0-based coordinates
# -ibam BAM file as input for coverage
#-g id the number of chromosome file and the length of each one (tab delimited)
#> your output file
bedtools genomecov -dz -ibam WPT391-1-6v2.sorted.bam \
-g /home/stuparr/michnoj0/Reference_Files/Soybean_chromosome.txt \
> WPT391-1-6v2_per_base_coverage.txt


########use Meesh R script to calculate mean and standard  deviation
module load R
Rscript Mean_and_SD_freqtable.R WPT391-1-6v2_per_base_coverage.txt > mean2SDWPT391-1-6v2.txt


#generate bcf file
#-c50 o reduce the effect of reads with excessive mismatches.
#-uf is the fasta of your referenct
#-bvcg to make a bcf file
#use bcf tools to turn it into a vcf
# use varFilter to -d for min read for SNP and -D for max read of SNP (mean+2SD)
samtools mpileup -c50 \
-uf /home/stuparr/shared/References/Gmax.a2.v1/assembly/Gmax_275_v2.0.fa \
WPT391-1-6v2.sorted.bam \
| bcftools view -bvcg - > WPT391-1-6v2_raw.bcf
bcftools view WPT391-1-6v2_raw.bcf > WPT391-1-6v2_raw.vcf
vcfutils.pl varFilter -d 3 -D 237 WPT391-1-6v2_raw.vcf > WPT391-1-6v2_Filter.vcf

#Analyze your data in R




















