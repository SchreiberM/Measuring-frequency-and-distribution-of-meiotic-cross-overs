#$/bin/bash

#Tools required: vcftools, bgzip, bcftools

#I/O
bcfFile=$1
sampleName=$2
bgName=$3
prefix=$4

## 1. bcftools to (1) keep only the SNPs
## (2) keep only those with allele balance (INFO/AB) between 0.25 to 0.75
## (3) Sample depth of mutant > 10 (FROMAT/DP)
## (4) Sample depth of Bowman > 25 (FORMAT/DP)
## (4) Maximum read depth < 70 (INFO/DP)
## (5) snpsift to to filter sites having GEN[BG]=Hom and GEN[mutant]=het
## Note. Need to keep non-biallelic sites too, because
## Bowman can have genotype different from the ref (ref is TT, Bowman is GG as 1/1), and a multiallelic site could also be a valid SNP (ex. the mutated is A/G as 1/2 )

## Then, carry out the filtering
vcftools --bcf $bcfFile --max-missing 1.0 --minQ 30 --minDP 10 --recode --recode-INFO-all --out $prefix.int.het
sed -i '1 i\##fileformat=VCFv4.2' $prefix.int.het.recode.vcf
bgzip -c  $prefix.int.het.recode.vcf > $prefix.int.het.recode.vcf.gz
tabix $prefix.int.het.recode.vcf.gz

bcftools view --types snps --include '(INFO/AB > 0.3 & INFO/AB < 0.7) & (INFO/AO > 2 ) & (INFO/DP < 100)' $prefix.int.het.recode.vcf.gz | java -jar $snpsift filter "isHet(GEN[${sampleName}]) & isHom(GEN[${bgName}])" > $prefix.snp.het.vcf


## 2. Using the vcf utils vcf-to-tab to convert the vcf to genotype tab format
## Thus the comparison between background and mutant can be done easily
## Keeping those where GEN[BG]=C/C and GEN[MUT]=C/T | T/C
## and                 GEN[BG]=G/G and GEN[MUT]=G/A | A/G
## Assuming BG genotype is on col-4 and MUT genotype is on col-5
cat $prefix.snp.het.vcf | vcf-to-tab > $prefix.snp.het.tab

## This step generates the $prefix.hetFiltered.tab for checking the genotype correctness
head -n 1 $prefix.snp.het.tab > $prefix.hetFiltered.tab
# Check which field (4 or 5) does the sample name locates
fieldName=$(awk '{print $4}' $prefix.hetFiltered.tab)
if   [ "$fieldName" = "$bgName" ];then
	awk '(($4=="G/G" && $5=="G/A") || ($4=="G/G" && $5=="A/G")) || (($4=="C/C" && $5=="C/T") || ($4=="C/C" && $5=="T/C"))' $prefix.snp.het.tab >> $prefix.hetFiltered.tab
elif [ "$fieldName" = "$sampleName" ]; then
	awk '(($5=="G/G" && $4=="G/A") || ($5=="G/G" && $4=="A/G")) || (($5=="C/C" && $4=="C/T") || ($5=="C/C" && $4=="T/C"))' $prefix.snp.het.tab >> $prefix.hetFiltered.tab
fi

## This step generates the .lst file for (1) plot snp density (2) later used for filtering by position (--positions in vcftools)
if   [ "$fieldName" = "$bgName" ];then
	awk -v bgName=$bgName -v sampleName=$sampleName -v OFS="\t" '{print $1,$2,$2,bgName":"$4","sampleName":"$5}' $prefix.hetFiltered.tab | grep -v "#"  > $prefix.hetFiltered.lst
elif [ "$fieldName" = "$sampleName" ]; then
	awk -v bgName=$bgName -v sampleName=$sampleName -v OFS="\t" '{print $1,$2,$2,bgName":"$5","sampleName":"$4}' $prefix.hetFiltered.tab | grep -v "#"  > $prefix.hetFiltered.lst
fi

## 2.5 Remove the introgressed regions for BW230.des10
if [ "$sampleName" = "BW230.des10.M2" ]; then
	## Rename the original (i.e., containing the introgressed region) lst to the original.lst file
	cp $prefix.hetFiltered.lst $prefix.hetFiltered.original.lst
	## Use bedtools subtract tool to remove the introgressed regions
	bedtools subtract -a $prefix.hetFiltered.original.lst -b introgressedRegions.bed > $prefix.hetFiltered.lst
	## Rename the original tab file to original.tab file
	mv $prefix.hetFiltered.tab $prefix.hetFiltered.original.tab
fi

## 3. Count how many snps there are after filtering
echo
echo ==========STATS OF $sampleName==========
echo Total number of SNPs - $(tail -n +2 $prefix.hetFiltered.lst | wc -l)
echo chr1H - $(grep "chr1H" $prefix.hetFiltered.lst | wc -l)
echo chr2H - $(grep "chr2H" $prefix.hetFiltered.lst | wc -l)
echo chr3H - $(grep "chr3H" $prefix.hetFiltered.lst | wc -l)
echo chr4H - $(grep "chr4H" $prefix.hetFiltered.lst | wc -l)
echo chr5H - $(grep "chr5H" $prefix.hetFiltered.lst | wc -l)
echo chr6H - $(grep "chr6H" $prefix.hetFiltered.lst | wc -l)
echo chr7H - $(grep "chr7H" $prefix.hetFiltered.lst | wc -l)
echo chrUn - $(grep "chrUn" $prefix.hetFiltered.lst | wc -l)
echo ========================================
echo

## 4. Extract to vcf file and compute stats
bcftools view --regions-file $prefix.hetFiltered.lst --output-type b --output-file $prefix.hetFiltered.bcf $bcfFile
bcftools index $prefix.hetFiltered.bcf
bcftools stats $prefix.hetFiltered.bcf > $prefix.hetFiltered.stats

