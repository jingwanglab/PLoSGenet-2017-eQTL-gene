#! /bin/bash -l

#SBATCH -A b2012141
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o eGene.core.gene_list.out
#SBATCH -e eGene.core.gene_list.err
#SBATCH -J eGene.core.gene_list.job
#SBATCH -t 2-00:00:00

module load bioinfo-tools

summary_perl="/proj/b2011141/pipeline/perl/eQTLs/eqtls.thetas.get_gene.pl"
attach_gene_perl="/proj/b2011141/pipeline/perl/eQTLs/eqtls.thetas.attach_gene.pl"
 

InputDir="/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas"
InputDir_fold="/proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/0_4_fold/new"

#for file in $InputDir/*pestPG
#do
#perl $summary_perl $file
#done

#for file in $InputDir/*gene.txt
#do
#sed 's/#//g' $file > $InputDir/temp && mv $InputDir/temp $file
#done

#for file in $InputDir/*thetas.txt
#do
#name=${file%%thetas.txt}
#perl $attach_gene_perl $file ${name}thetas.gz.gene.txt
#done

###############0-fold, 4-fold###################
#for file in $InputDir_fold/*pestPG
#do
#perl $summary_perl $file
#done

for file in $InputDir_fold/*gene.txt
do
sed 's/#//g' $file > $InputDir_fold/temp && mv $InputDir_fold/temp $file
done

for file in $InputDir_fold/*thetas.txt
do
name=${file%%thetas.txt}
perl $attach_gene_perl $file ${name}thetas.gz.gene.txt
done


