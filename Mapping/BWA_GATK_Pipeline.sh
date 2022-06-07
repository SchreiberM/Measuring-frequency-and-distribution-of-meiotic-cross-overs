#!/bin/bash

#$ -V

##########################################################################################################
#This is a pipeline script for aligning PE Illumina reads from a single sample 
#to a FASTA reference sequence and then calling variants on it using the GATK.
#It follows the GATK best practice in broad terms and assumes there is no existing
#truth data available in terms of known SNPs (as is the case for most non-model organisms). It
#creates its own truth data in an initial run of the HaplotypeCaller and uses these as the basis for 
#the base quality score recalibration. 

#Tools required: BWA, sambamba, vcflib, GATK3

##########################################################################################################
#variables -- mostly from command line
##########################################################################################################
#check for the correct number of args
if [ $# -ne 6 ]
then
    echo "Error in $0 - Invalid Argument Count"
    echo "Syntax: $0  <R1 FASTQ file, gzipped> <R2 FASTQ file, gzipped> <sample name> <full path to reference sequence, FASTA> <min. alignment score> <numThreads>"
    exit
fi

#R1 FASTQ file, gzipped
R1=$1
#R2 FASTQ file, gzipped
R2=$2
#sample name -- can be different from file name
sample=$3
#full path to reference sequence, FASTA
refseq=$4
#min. alignment score
numberMismatches=$5
#number of threads
numThreads=$6


###########################################################################################################
##map the reads with BWA-MEM
###########################################################################################################

echo -e "\n========================================="
echo "mapping with BWA"
echo "start time `/bin/date`"
echo -e "========================================="

bwa mem \
$refseq \
<(pigz -cd $R1) \
<(pigz -cd $R2) \
-t $numThreads \
-R "@RG\tID:$sample\tSM:$sample\tPL:Illumina" 2> $sample.bwa.log \
> $sample.sam

echo -e "\n========================================="
echo "convert SAM to BAM"
echo "start time `/bin/date`"
echo -e "========================================="

sambamba view \
--sam-input \
--format=bam \
-l 0 \
-t $numThreads \
--filter="[NM] <= $numberMismatches" \
-o $sample.unsorted.bam \
$sample.sam

rm $sample.sam

##########################################################################################################
#Sort BAM file
##########################################################################################################
echo -e "\n========================================="
echo "sorting BAM file"
echo "start time `/bin/date`"
echo -e "========================================="

#echo "sambamba sort -t $numThreads -o $sample.sorted.bam $sample.unsorted.bam"
sambamba sort -t $numThreads -o $sample.sorted.bam $sample.unsorted.bam

#echo "samtools sort $sample.unsorted.bam  $sample.sorted"  # note: samtools 1.2
samtools sort $sample.unsorted.bam -o $sample.sorted 
rm $sample.unsorted.bam

##########################################################################################################
#duplicate removal
##########################################################################################################
echo -e "\n========================================="
echo "removing duplicates"
echo "start time `/bin/date`"
echo -e "========================================="

sambamba markdup --overflow-list-size 600000 -r -t $numThreads $sample.sorted $sample.rmduped.bam 
rm $sample.sorted.bam

##########################################################################################################
#indel realignment
##########################################################################################################
echo -e "\n========================================="
echo "Indel realignment step 1 - generating target interval list"
echo "start time `/bin/date`"
echo -e "========================================="

gatk3 -T RealignerTargetCreator -R $refseq -I $sample.rmduped.bam -o $sample.target_intervals.list -nt $numThreads

echo -e "\n========================================="
echo "Indel realignment step 2 - realigning"
echo "start time `/bin/date`"
echo -e "========================================="

gatk3 -T IndelRealigner -R $refseq -I $sample.rmduped.bam -targetIntervals $sample.target_intervals.list -o $sample.realigned.bam

#index the newly created realigned BAM file
sambamba index -t $numThreads $sample.realigned.bam $sample.realigned.bam.bai

rm $sample.rmduped.bam
rm $sample.rmduped.bam.bai

##########################################################################################################
#initial run of the HaplotypeCaller -- this generates the first set of variants we will use as truth data for the base quality score recalibration
##########################################################################################################
echo -e "\n========================================="
echo "Running the HaplotypeCaller - first run"
echo "start time `/bin/date`"
echo -e "========================================="

gatk3 -T HaplotypeCaller -R $refseq -I $sample.realigned.bam -o $sample.initial.variants.vcf -dontUseSoftClippedBases -nct $numThreads

##########################################################################################################
#filter the initial variants to remove poor quality calls with a QUAL of
#less than 20 - this will become our truth set for the BQSR step
##########################################################################################################
echo "filtering variants for >QUAL20"

vcffilter -f "QUAL > 20" $sample.initial.variants.vcf > $sample.initial.variants.filteredQ20.vcf

##########################################################################################################
#base quality score recalibration
##########################################################################################################
#produce the recalibration table
echo -e "\n========================================="
echo "producing the recalibration table"
echo "start time `/bin/date`"
echo -e "========================================="

gatk3 -T BaseRecalibrator -R $refseq -I $sample.realigned.bam -knownSites $sample.initial.variants.filteredQ20.vcf -o $sample.recal.table -nct $numThreads

#having produced the table, now recalibrate the data
echo -e "\n========================================="
echo "recalibrating the data"
echo "start time `/bin/date`"
echo -e "========================================="

gatk3 -T PrintReads -R $refseq -I $sample.realigned.bam -BQSR $sample.recal.table -o $sample.recalibrated.bam -nct $numThreads

rm $sample.realigned.bai
rm $sample.realigned.bam
rm $sample.realigned.bam.bai

##########################################################################################################
#done
##########################################################################################################
echo -e "\n========================================="
echo "PIPELINE COMPLETE"
echo "END TIME `/bin/date`"
echo -e "========================================="

