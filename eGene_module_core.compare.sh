#! /bin/bash -l

#SBATCH -A b2011141
#SBATCH -p core
#SBATCH -n 2
#SBATCH -o eGene_vs_noneGene.compare.out
#SBATCH -e eGene_vs_noneGene.compare.out
#SBATCH -J eGene_vs_noneGene.compare.out
#SBATCH -t 2-00:00:00

module load bioinfo-tools
module load BEDTools


###This script is used to compare the paramters used to estimate the extent of selection for eGenes and for those non-Genes, in addition, also compare the module core genes or non-core genes in each eGene or non-eGene group

###The parameters we want to test are: theta, TajimaD, theta_0fold/theta_4fold, DFE, d_0fold/d_4fold

###1. Extract the 4-fold and 0-fold sites from the .gff annotation file

#zcat /proj/b2011141/nobackup/gff/tremula/v1.1/gff3/Potra01-gene-representative.gff3.gz | grep -v "#" | grep "CDS" |cut -f 1,4,5,7,9 |cut -d ";" -f 1 |cut -d "." -f 1 |sed 's/ID=//g' > /proj/b2011141/pipeline/python/tremula/tremula.cds.longest.txt

#python AnnotateRef.py /proj/b2010014/GenomePaper/Populus-tremula/v1.0/v1.0/fasta/Potra01-genome.fa tremula.cds.longest.txt > tremula.cds.longest.annotation.txt

###2. Esimate the mean theta and Tajima's D for genes in the four groups

###transfer thetas.gz file to bed file in order to use BEDTools

#perl /proj/b2011141/pipeline/perl/thetas_to_bed.pl /proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/ANGSD/SFS/all/total/thetas/SwAsp_94.total.thetas.gz

###group 1: Genes that are both eGenes and core genes in the module

#awk '$2=="TRUE"' /proj/b2011141/nobackup/eQTL_paper/gene_list/swasp_gene_table.txt |awk '$4=="TRUE"' | awk '$5=="TRUE"' |cut -f 1 > /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.txt
#zcat /proj/b2011141/nobackup/gff/tremula/v1.1/gff3/Potra01-gene-representative.gff3.gz |grep -v "#" | grep "gene" | grep -f /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.txt > /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.gff.txt && rm /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.txt

###group 2: Genes that are both eGenes and non-core genes in the module

#awk '$2=="TRUE"' /proj/b2011141/nobackup/eQTL_paper/gene_list/swasp_gene_table.txt |awk '$4=="FALSE"' | awk '$5=="TRUE"' |cut -f 1 > /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.txt
#zcat /proj/b2011141/nobackup/gff/tremula/v1.1/gff3/Potra01-gene-representative.gff3.gz |grep -v "#" | grep "gene" | grep -f /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.txt > /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.gff.txt && rm /proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.txt

###group 3: Genes that are not eGenes but core genes in the module

#awk '$2=="TRUE"' /proj/b2011141/nobackup/eQTL_paper/gene_list/swasp_gene_table.txt |awk '$4=="TRUE"' | awk '$5=="FALSE"' |cut -f 1 > /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.txt
#zcat /proj/b2011141/nobackup/gff/tremula/v1.1/gff3/Potra01-gene-representative.gff3.gz |grep -v "#" | grep "gene" | grep -f /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.txt > /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.gff.txt && rm /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.txt

###group4: Genes that are not eGenes and not core genes in the module

#awk '$2=="TRUE"' /proj/b2011141/nobackup/eQTL_paper/gene_list/swasp_gene_table.txt |awk '$4=="FALSE"' | awk '$5=="FALSE"' |cut -f 1 > /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.txt
#zcat /proj/b2011141/nobackup/gff/tremula/v1.1/gff3/Potra01-gene-representative.gff3.gz |grep -v "#" | grep "gene" | grep -f /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.txt > /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.gff.txt && rm /proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.txt


OutDir="/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas"

angsd="/proj/b2011141/tools/angsd0.902/angsd/angsd"
realSFS="/proj/b2011141/tools/angsd0.902/angsd/misc/realSFS"
thetaStat="/proj/b2011141/tools/angsd0.902/angsd/misc/thetaStat"
nChrom=94

thetas_bed="/proj/b2011141/nobackup/PaperIII-local_adaptation/asp201_94/ANGSD/SFS/all/total/thetas/SwAsp_94.total.bed"

####eGene_core
#eGene_core="/proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.gff.txt"
#n_eGene_core=$(cat $eGene_core | wc -l)


#bedtools intersect -a $thetas_bed -b $eGene_core -wb > $OutDir/eGene_core.gff.thetas.bed

#for i in `seq 1 $n_eGene_core`; do 
#gene=$(head -n $i $eGene_core | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir/eGene_core.gff.thetas.bed | cut -f 1,3,4,5,6,7,8 > $OutDir/${gene}.thetas
#gzip $OutDir/${gene}.thetas
#$thetaStat make_bed $OutDir/${gene}.thetas.gz
#$thetaStat do_stat $OutDir/${gene}.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir/eGene_core.thetas.gz.pestPG
#cat $OutDir/${gene}.thetas.gz.pestPG >> $OutDir/eGene_core.thetas.gz.pestPG && rm $OutDir/${gene}.thetas.gz $OutDir/${gene}.thetas.gz.bin $OutDir/${gene}.thetas.gz.idx $OutDir/${gene}.thetas.gz.pestPG
#done

###eGene_noncore
#eGene_noncore="/proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.gff.txt"
#n_eGene_noncore=$(cat $eGene_noncore | wc -l)

#bedtools intersect -a $thetas_bed -b $eGene_noncore -wb > $OutDir/eGene_noncore.gff.thetas.bed
#
#for i in `seq 1 $n_eGene_noncore`; do
#gene=$(head -n $i $eGene_noncore | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir/eGene_noncore.gff.thetas.bed | cut -f 1,3,4,5,6,7,8 > $OutDir/${gene}.thetas
#gzip $OutDir/${gene}.thetas
#$thetaStat make_bed $OutDir/${gene}.thetas.gz
#$thetaStat do_stat $OutDir/${gene}.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir/eGene_noncore.thetas.gz.pestPG
#cat $OutDir/${gene}.thetas.gz.pestPG >> $OutDir/eGene_noncore.thetas.gz.pestPG && rm $OutDir/${gene}.thetas.gz $OutDir/${gene}.thetas.gz.bin $OutDir/${gene}.thetas.gz.idx $OutDir/${gene}.thetas.gz.pestPG
#done

###non_eGene_core
#non_eGene_core="/proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.gff.txt"
#n_non_eGene_core=$(cat $non_eGene_core | wc -l)

#bedtools intersect -a $thetas_bed -b $non_eGene_core -wb > $OutDir/non_eGene_core.gff.thetas.bed

#for i in `seq 1 $n_non_eGene_core`; do
#gene=$(head -n $i $non_eGene_core | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir/non_eGene_core.gff.thetas.bed | cut -f 1,3,4,5,6,7,8 > $OutDir/${gene}.thetas
#gzip $OutDir/${gene}.thetas
#$thetaStat make_bed $OutDir/${gene}.thetas.gz
#$thetaStat do_stat $OutDir/${gene}.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir/non_eGene_core.thetas.gz.pestPG
#cat $OutDir/${gene}.thetas.gz.pestPG >> $OutDir/non_eGene_core.thetas.gz.pestPG && rm $OutDir/${gene}.thetas.gz $OutDir/${gene}.thetas.gz.bin $OutDir/${gene}.thetas.gz.idx $OutDir/${gene}.thetas.gz.pestPG
#done

###non_eGene_noncore
#non_eGene_noncore="/proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.gff.txt"
#n_non_eGene_noncore=$(cat $non_eGene_noncore | wc -l)


#bedtools intersect -a $thetas_bed -b $non_eGene_noncore -wb > $OutDir/non_eGene_noncore.gff.thetas.bed

#for i in `seq 1 $n_non_eGene_noncore`; do
#gene=$(head -n $i $non_eGene_noncore | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir/non_eGene_noncore.gff.thetas.bed | cut -f 1,3,4,5,6,7,8 > $OutDir/${gene}.thetas
#gzip $OutDir/${gene}.thetas
#$thetaStat make_bed $OutDir/${gene}.thetas.gz
#$thetaStat do_stat $OutDir/${gene}.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir/non_eGene_noncore.thetas.gz.pestPG
#cat $OutDir/${gene}.thetas.gz.pestPG >> $OutDir/non_eGene_noncore.thetas.gz.pestPG && rm $OutDir/${gene}.thetas.gz $OutDir/${gene}.thetas.gz.bin $OutDir/${gene}.thetas.gz.idx $OutDir/${gene}.thetas.gz.pestPG
#done


#####3. Esimating the theta and Tajima's D separately for zero-fold and four-fold sites

OutDir_fold="/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/0_4_fold/new"

###use perl script to transfer the .txt file to .bed file
#perl /proj/b2011141/pipeline/perl/0_4_fold_to_bed.pl /proj/b2011141/pipeline/python/tremula/0_fold.txt 
#perl /proj/b2011141/pipeline/perl/0_4_fold_to_bed.pl /proj/b2011141/pipeline/python/tremula/4_fold.txt 

zero_fold_bed="/proj/b2011141/pipeline/python/tremula/0_fold.bed"
four_fold_bed="/proj/b2011141/pipeline/python/tremula/4_fold.bed"


####3.1eGene_core
#bedtools intersect -a $zero_fold_bed -b $OutDir/eGene_core.gff.thetas.bed -wb > $OutDir_fold/eGene_core.zero_fold.thetas.bed
#bedtools intersect -a $four_fold_bed -b $OutDir/eGene_core.gff.thetas.bed -wb > $OutDir_fold/eGene_core.four_fold.thetas.bed

eGene_core="/proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_core.gff.txt"
n_eGene_core=$(cat $eGene_core | wc -l)


#for i in `seq 1 $n_eGene_core`; do
#gene=$(head -n $i $eGene_core | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir_fold/eGene_core.zero_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.zero_fold.thetas
#grep "$gene" $OutDir_fold/eGene_core.four_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.four_fold.thetas
#gzip $OutDir_fold/${gene}.zero_fold.thetas
#gzip $OutDir_fold/${gene}.four_fold.thetas
#$thetaStat make_bed $OutDir_fold/${gene}.zero_fold.thetas.gz
#$thetaStat make_bed $OutDir_fold/${gene}.four_fold.thetas.gz
#$thetaStat do_stat $OutDir_fold/${gene}.zero_fold.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir_fold/${gene}.four_fold.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir_fold/eGene_core.zero_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG  >> $OutDir_fold/eGene_core.zero_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.zero_fold.thetas.gz $OutDir_fold/${gene}.zero_fold.thetas.gz.bin $OutDir_fold/${gene}.zero_fold.thetas.gz.idx $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG
#echo -e "#$gene" >> $OutDir_fold/eGene_core.four_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG  >> $OutDir_fold/eGene_core.four_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.four_fold.thetas.gz $OutDir_fold/${gene}.four_fold.thetas.gz.bin $OutDir_fold/${gene}.four_fold.thetas.gz.idx $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG
#done


###3.2 eGene_noncore
eGene_noncore="/proj/b2011141/nobackup/eQTL_paper/gene_list/eGene_noncore.gff.txt"
n_eGene_noncore=$(cat $eGene_noncore | wc -l)

#bedtools intersect -a $zero_fold_bed -b $OutDir/eGene_noncore.gff.thetas.bed -wb > $OutDir_fold/eGene_noncore.zero_fold.thetas.bed
#bedtools intersect -a $four_fold_bed -b $OutDir/eGene_noncore.gff.thetas.bed -wb > $OutDir_fold/eGene_noncore.four_fold.thetas.bed

#for i in `seq 1 $n_eGene_noncore`; do
#gene=$(head -n $i $eGene_noncore | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir_fold/eGene_noncore.zero_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.zero_fold.thetas
#grep "$gene" $OutDir_fold/eGene_noncore.four_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.four_fold.thetas
#gzip $OutDir_fold/${gene}.zero_fold.thetas
#gzip $OutDir_fold/${gene}.four_fold.thetas
#$thetaStat make_bed $OutDir_fold/${gene}.zero_fold.thetas.gz
#$thetaStat make_bed $OutDir_fold/${gene}.four_fold.thetas.gz
#$thetaStat do_stat $OutDir_fold/${gene}.zero_fold.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir_fold/${gene}.four_fold.thetas.gz -nChr $nChrom
#echo -e "#$gene" >> $OutDir_fold/eGene_noncore.zero_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG >> $OutDir_fold/eGene_noncore.zero_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.zero_fold.thetas.gz $OutDir_fold/${gene}.zero_fold.thetas.gz.bin $OutDir_fold/${gene}.zero_fold.thetas.gz.idx $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG
#echo -e "#$gene" >> $OutDir_fold/eGene_noncore.four_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG >> $OutDir_fold/eGene_noncore.four_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.four_fold.thetas.gz $OutDir_fold/${gene}.four_fold.thetas.gz.bin $OutDir_fold/${gene}.four_fold.thetas.gz.idx $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG
#done

###3.3 non_eGene_core

#bedtools intersect -a $zero_fold_bed -b $OutDir/non_eGene_core.gff.thetas.bed -wb > $OutDir_fold/non_eGene_core.zero_fold.thetas.bed
#bedtools intersect -a $four_fold_bed -b $OutDir/non_eGene_core.gff.thetas.bed -wb > $OutDir_fold/non_eGene_core.four_fold.thetas.bed

non_eGene_core="/proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_core.gff.txt"
n_non_eGene_core=$(cat $non_eGene_core | wc -l)


#for i in `seq 1 $n_non_eGene_core`; do
#gene=$(head -n $i $non_eGene_core | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir_fold/non_eGene_core.zero_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.zero_fold.thetas
#grep "$gene" $OutDir_fold/non_eGene_core.four_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.four_fold.thetas
#gzip $OutDir_fold/${gene}.zero_fold.thetas
#gzip $OutDir_fold/${gene}.four_fold.thetas
#
#$thetaStat make_bed $OutDir_fold/${gene}.zero_fold.thetas.gz
#$thetaStat make_bed $OutDir_fold/${gene}.four_fold.thetas.gz
#$thetaStat do_stat $OutDir_fold/${gene}.zero_fold.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir_fold/${gene}.four_fold.thetas.gz -nChr $nChrom
#
#echo -e "#$gene" >> $OutDir_fold/non_eGene_core.zero_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG >> $OutDir_fold/non_eGene_core.zero_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.zero_fold.thetas.gz $OutDir_fold/${gene}.zero_fold.thetas.gz.bin $OutDir_fold/${gene}.zero_fold.thetas.gz.idx $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG
#
#echo -e "#$gene" >> $OutDir_fold/non_eGene_core.four_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG >> $OutDir_fold/non_eGene_core.four_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.four_fold.thetas.gz $OutDir_fold/${gene}.four_fold.thetas.gz.bin $OutDir_fold/${gene}.four_fold.thetas.gz.idx $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG
#done


###3.4 non_eGene_noncore

#bedtools intersect -a $zero_fold_bed -b $OutDir/non_eGene_noncore.gff.thetas.bed -wb > $OutDir_fold/non_eGene_noncore.zero_fold.thetas.bed
#bedtools intersect -a $four_fold_bed -b $OutDir/non_eGene_noncore.gff.thetas.bed -wb > $OutDir_fold/non_eGene_noncore.four_fold.thetas.bed

non_eGene_noncore="/proj/b2011141/nobackup/eQTL_paper/gene_list/non_eGene_noncore.gff.txt"
n_non_eGene_noncore=$(cat $non_eGene_noncore | wc -l)

#for i in `seq 1 $n_non_eGene_noncore`; do
#
#gene=$(head -n $i $non_eGene_noncore | tail -n 1|cut -f 9 | cut -f 1 -d ";" |sed 's/ID=//g')
#grep "$gene" $OutDir_fold/non_eGene_noncore.zero_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.zero_fold.thetas
#grep "$gene" $OutDir_fold/non_eGene_noncore.four_fold.thetas.bed | cut -f 1,3,7,8,9,10,11 > $OutDir_fold/${gene}.four_fold.thetas
#
#gzip $OutDir_fold/${gene}.zero_fold.thetas
#gzip $OutDir_fold/${gene}.four_fold.thetas
#
#$thetaStat make_bed $OutDir_fold/${gene}.zero_fold.thetas.gz
#$thetaStat make_bed $OutDir_fold/${gene}.four_fold.thetas.gz
#$thetaStat do_stat $OutDir_fold/${gene}.zero_fold.thetas.gz -nChr $nChrom
#$thetaStat do_stat $OutDir_fold/${gene}.four_fold.thetas.gz -nChr $nChrom
#
#echo -e "#$gene" >> $OutDir_fold/non_eGene_noncore.zero_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG >> $OutDir_fold/non_eGene_noncore.zero_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.zero_fold.thetas.gz $OutDir_fold/${gene}.zero_fold.thetas.gz.bin $OutDir_fold/${gene}.zero_fold.thetas.gz.idx $OutDir_fold/${gene}.zero_fold.thetas.gz.pestPG
#
#echo -e "#$gene" >> $OutDir_fold/non_eGene_noncore.four_fold.thetas.gz.pestPG
#cat $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG >> $OutDir_fold/non_eGene_noncore.four_fold.thetas.gz.pestPG && rm $OutDir_fold/${gene}.four_fold.thetas.gz $OutDir_fold/${gene}.four_fold.thetas.gz.bin $OutDir_fold/${gene}.four_fold.thetas.gz.idx $OutDir_fold/${gene}.four_fold.thetas.gz.pestPG
#
#done


###remove the gene that do not have thetas data
#grep -v '^Potra' $OutDir_fold/eGene_core.four_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/eGene_core.four_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/eGene_noncore.four_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/eGene_noncore.four_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/non_eGene_core.four_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/non_eGene_core.four_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/non_eGene_noncore.four_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/non_eGene_noncore.four_fold.thetas.gz.pestPG

#grep -v '^Potra' $OutDir_fold/eGene_core.zero_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/eGene_core.zero_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/eGene_noncore.zero_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/eGene_noncore.zero_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/non_eGene_core.zero_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/non_eGene_core.zero_fold.thetas.gz.pestPG
#grep -v '^Potra' $OutDir_fold/non_eGene_noncore.zero_fold.thetas.gz.pestPG > $OutDir_fold/temp && mv $OutDir_fold/temp $OutDir_fold/non_eGene_noncore.zero_fold.thetas.gz.pestPG


######4. R statistics estimates

###4.1.1 mean theta,pi, Tajima's D separately for the four groups
eQTL_thetas="/proj/b2011141/pipeline/R/eQTLs/eQTL.thetas_table.R"

#Rscript $eQTL_thetas $OutDir eGene_core
#Rscript $eQTL_thetas $OutDir eGene_noncore
#Rscript $eQTL_thetas $OutDir non_eGene_core
#Rscript $eQTL_thetas $OutDir non_eGene_noncore
####4.1.2 make the plot
#Rscript /proj/b2011141/pipeline/R/eQTLs/eQTL.thetas_plot.R


###4.2.1 theta 0_fold/4_fold separately for each of the four groups

eQTL_thetas_fold="/proj/b2011141/pipeline/R/eQTLs/eQTL.thetas_table.0_4fold.R"
Rscript $eQTL_thetas_fold $OutDir_fold eGene_core
Rscript $eQTL_thetas_fold $OutDir_fold eGene_noncore
Rscript $eQTL_thetas_fold $OutDir_fold non_eGene_core
Rscript $eQTL_thetas_fold $OutDir_fold non_eGene_noncore




