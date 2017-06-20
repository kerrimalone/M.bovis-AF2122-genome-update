#############################################
#    Unmapped read assembly: Identification #
#		        of RD900                    #
# --- Linux bioinformatics workflow ---     #
#############################################
# Author of current version (1.0.0): Malone, K.M

##################################
# Data download and md5sum check #
##################################

#Download data using appropriate protocol, username and password
cd $HOME/storage/DNAseq/bovis_AF2122
mkdir raw_data
cd !$
nohup wget -r ftp<server> --ftp_user='<username>' --ftp-password='<password>' &

#Check the md5sum files to insure successful transfer and file integrity
mkdir -p $HOME/storage/DNAseq/bovis_AF2122/md5check
cd !$
for file in 'find $HOME/storage/DNAseq/bovis_AF2122 -name *fastq.gz'; \
do echo "md5sum $file >> md5sum_check.txt"; \
done; 

#Split the script and run
split -d -l 92 md5sum.sh md5sum.sh.
for script in `ls md5sum.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

#Check output and proceed

##########################################
# FastQC quality check of raw FASTQ files #
###########################################

# Required software is FastQC v0.11.5, consult manual/tutorial
# for details: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Create and enter the quality check output directory:
mkdir -p $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter
cd !$

# Run FastQC in one file to see if it's working well:
fastqc -o $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter \
--noextract --nogroup -t 2 \
$HOME/storage/DNAseq/bovis_AF2122/raw_data/2122-1_S9_L001_R1_001.fastq.gz

#Check output and proceed

# Create a bash script to perform FastQC quality check on all fastq.gz files:
for file in `find $HOME/storage/DNAseq/bovis_AF2122/raw_data/ \
-name *fastq.gz`; do echo "fastqc --noextract --nogroup -t 1 \
-o $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter $file" \
>> fastqc.sh; done;

#Split the script and run
split -d -l 70 fastqc.sh fastqc.sh.
for script in `ls fastqc.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

# Check if all of the files were processed:
for file in `ls fastqc.sh.0*.nohup`; \
do more $file | grep "Failed to process file" >> failed_fastqc.txt
done

# Deleted all the HTML files:
rm -r *.html

# Check all output from FastQC:
mkdir $HOME/scratch/PPDbRNAseqTimeCourse/quality_check/pre-filtering/tmp

for file in `ls *_fastqc.zip`; do unzip \
$file -d $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter/tmp; \
done;

for file in \
`find $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter/tmp \
-name summary.txt`; do more $file >> reports_pre-filtering.txt; done

for file in \
`find $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter/tmp \
-name fastqc_data.txt`; do head -n 10 $file >> basic_stats_pre-filtering.txt; \
done

# Check outputs, remove temporary folder and the associated files:
rm -rf $HOME/scratch/DNAseq/bovis_AF2122/QC/pre_filter/tmp

##########################################
#          Read trimming    			 #
##########################################
#Required software is Sickle,
#see https://github.com/najoshi/sickle for more information and download.


#Trim reads based on quality score from fastQC (per base quality sequence graph) with Sickle. Using a sliding window quality cut-off score of 30.
for file in `find $HOME/storage/DNAseq/bovis_AF2122/raw_data/ \
-name *fastq.gz`; 
do file2=`echo $file | perl -p -e 's/(fastq.gz)/_q30.fastq.gz/'`; \
do file3=`echo $file | perl -p -e 's/(fastq.gz)/_singles.fastq.gz/'`; \
echo "sickle pe -c $file -t sanger -m $file2 -s $file3 -q30 -l15" \
>> trimming.sh; \
done;

#Split the script and run
split -d -l 70 trimming.sh trimming.sh.
for script in `ls trimming.sh.*`
do
chmod 755 $script
nohup ./$script > ${script}.nohup &
done

#Re-run fastQC on the trimmed files and ensure trimming was successful


########################
# Mapping to BX28333.2 #
########################
#Mapping Illumina reads to the M. bovis AF2122/97 (BX28333.2) genome sequence
#Software required is Stampy v1.0.20,
#more details at http://www.well.ox.ac.uk/~gerton/README.txt
#Merge replicates but keep read mates separated
cd $HOME/storage/DNAseq/bovis_AF2122/raw_data/
cat 2122-1_S9_L001_R1_001*q30.fastq.gz 2122-2_S12_L001_R1_001*q30.fastq.gz > 2122_R1.fastq.gz
cat 2122-1_S9_L001_R2_001*q30.fastq.gz 2122-2_S12_L001_R2_001*q30.fastq.gz > 2122_R2.fastq.gz

mkdir $HOME/scratch/DNAseq/bovis_AF2122/mapping
cd !$

#Build Stampy index and hash files based on de novo assembly
stampy -G BX28333.2 $HOME/scratch/DNAseq/genome_files/BX28333.2.fasta
stampy -g BX28333.2 -H BX28333.2

#Trimmed Illumina data used. For more information see de_novo_assembly.sh
nohup stampy -g $HOME/scratch/DNAseq/genome_files/BX28333.2 -h $HOME/scratch/DNAseq/genome_files/BX28333.2 \
-t10 -M $HOME/storage/DNAseq/bovis_AF2122/raw_data/2122_R1.fastq.gz,$HOME/storage/DNAseq/bovis_AF2122/raw_data/2122_R2.fastq.gz \
-o BX28333.2_illumina_mapped.sam > BX28333.2_illumina_mapped_summary.txt &

#Convert output sam to bam file
nohup samtools faidx $HOME/scratch/DNAseq/genome_files/BX28333.2.fasta &
nohup samtools import $HOME/scratch/DNAseq/genome_files/BX28333.2.fasta.fai BX28333.2_illumina_mapped.sam BX28333.2_illumina_mapped.bam &

#Get unmapped reads from output file
nohup samtools view -b -f 4 BX28333.2_illumina_mapped.bam > BX28333.2_illumina_unmapped.bam &

##############################
# Assembly of unmapped reads #
##############################

#Edit .config file for appropriate parameters
#maximal read length
max_rd_len=50 
[LIB]
avg_ins=50
#illumina data, forward-reverse
reverse_seq=1
asm_flags=3
rank=1
pair_num_cutoff=5
map_len=35
#path to fasta file containing both mates. You can input bams or fastqs but these need to be read mate separated files
p=$HOME/scratch/DNAseq/bovis_AF2122/mapping/BX28333.2_illumina_unmapped.bam

#Run SOAPdenovo-63mer
nohup ./SOAPdenovo-63mer all -s file.config -K 63 -R -p 10 -o BX28333.2_illumina_unmapped 1>ass.log 2>ass.err &

# Output file BX28333.2_illumina_unmapped.scafSeq contains assembled contigs.  
# Contigs were queried using BLAST and RD900 was identified.



