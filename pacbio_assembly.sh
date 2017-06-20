###########################################
#    De novo assembly of PacBio data      #
# --- Linux bioinformatics workflow ---   #
###########################################
# Author of current version (1.0.0): Malone, K.M

#Raw data was downloaded to $HOME/storage/DNAseq/Pacbio 
#and consisted of filtered subreads fastq files

#2 SMRT cells
#Min. polymerase read quality = 75
# Min. read length = 50 bp
# Min. subread length = 50 bp
# Mean subread length = 5,006 bp
# Subread N50 = 6, 339 bp
# Subread total number of bases = 539, 906, 876 bp


####################
# de novo assembly #
####################
#Required software is Canu,
#more details found here https://github.com/marbl/canu

#Merge filtered fastq files
cd $HOME/scratch/DNAseq/Pacbio
cat m160927_211843_42153_*subreads.fastq.gz > AF2122_filtered_subreads.fastq

mdkir $HOME/kmalone/scratch/DNAseq/Pacbio
cd !$

#Run Canu assembler using pre-optimised error rate = 0.013 (default)
./canu \
 -p AF2122 -d $HOME/scratch/DNAseq/Pacbio \
 genomeSize=4.4m \
 -pacbio-raw $HOME/storage/DNAseq/Pacbio/AF2122_filtered_subreads.fastq \
 gnuplotTested=true
 masterMemory 20 \
 masterThreads 10

######################
# Polishing assembly #
######################
#Map Illumina data to de novo assembly and polish

#Mapping Illumina reads to the de novo Pacbio assembly.
#Software required is Stampy v1.0.20,
#more details at http://www.well.ox.ac.uk/~gerton/README.txt

#Build Stampy index and hash files based on de novo assembly
stampy -G AF2122 $HOME/scratch/DNAseq/Pacbio/AF2122.contigs.fasta
stampy -g AF2122 -H AF2122

mkdir $HOME/scratch/DNAseq/Pacbio/illumina_mapping
cd!$
#Trimmed Illumina data used. For more information see de_novo_assembly.sh
nohup stampy -g ../AF2122 -h ../AF2122  -t10 -M $HOME/storage/DNAseq/bovis_AF2122/raw_data/2122_R1.fastq.gz,$HOME/storage/DNAseq/bovis_AF2122/raw_data/2122_R2.fastq.gz \
-o AF2122_illumina_mapped.sam > AF2122_illumina_mapped_summary.txt &


#Required software is Pilon v1.21,
#more details at https://github.com/broadinstitute/pilon/releases

#Create SAMtools index for Pacbio genome 
samtools faidx $HOME/scratch/DNAseq/Pacbio/AF2122.contigs.fasta

#Convert sam-to-bam
for file in `find $HOME/scratch/DNAseq/Pacbio/illumina_mapping \
-name *.sam`; \
do file2=`echo $file | perl -p -e 's/(.sam)/.bam/'`; \
echo "samtools import $HOME/scratch/DNAseq/Pacbio/*fasta.fai $file $file2" \
done;

#Sort bam
for file in `find $HOME/scratch/DNAseq/Pacbio/illumina_mapping \
-name *.bam`; \
do file2=`echo $file | perl -p -e 's/(.bam)/.sorted_bam/'`; \
echo "samtools sort $file $file2" \
done;

mkdir $HOME/scratch/DNAseq/Pacbio/pilon
cd !$
#Polishing the std out pacbio de novo assemblies with illumina data.
nohup java -Xmx16G -jar pilon-1.20.jar \
--genome ~$HOME/scratch/DNAseq/Pacbio/AF2122.contigs.fasta \
--bam ~$HOME/scratch/DNAseq/Pacbio/illumina_mapping/AF2122_illumina_mapped_sorted.bam  \
--output AF2122 --outdir $HOME/scratch/DNAseq/Pacbio/pilon --changes \
--vcf --tracks > pilon_summary.txt &


####################################
# Genome alignment using NUCmer    #
#        Pacbio assembly           #
#     versus AF2122/97(BX28333.2)  #
####################################
#Required software is MUMMER 3+,
#more details at http://mummer.sourceforge.net/manual/#nucmer

./nucmer --prefix=BX_AF2122_comparison $HOME/scratch/DNAseq/Pacbio/genome_files/BX248333.2.fasta $HOME/scratch/DNAseq/Pacbio/pilon/AF2122.fasta 
./show-snps -Clr BX_AF2122_comparison.delta > BX_AF2122_comparison.snps






