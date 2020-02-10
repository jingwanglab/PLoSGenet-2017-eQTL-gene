#! /bin/bash -l

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 3
#SBATCH -o eGene_gkaks.blat.out
#SBATCH -e eGene_gkaks.blat.err
#SBATCH -J eGene_gkaks.blat.job
#SBATCH -t 3-00:00:00

module load bioinfo-tools
module load blat
module load blast
module load paml
##first only select the genes in the Potri gff3 file
##grep -v "gene" Ptrichocarpa_210_v3.0.gene.gff3 |cut -f 1 -d ";" |sed 's/.v3.0.*//g' |sed 's/ID=//g' > /proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.gff

OutDir="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/gkaks"
gkaks="/proj/b2011141/tools/gKaKs/gKaKs_v1.3.pl"

###The difference between this and the other one is to let the blat to find the othologus sequences by itself
#perl $gkaks -query_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.cds.fa" -gff="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.gff" -hit_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potra.cds.fa" -spe=6 -blast2seq="bl2seq" -program="codeml" -blat="blat" -sfbh=T -gene_pair="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.Potra.txt" -kaks_file="$OutDir/potra_potri.kaks.txt" -detail="$OutDir/potra_potri.kaks.detail.txt"
perl $gkaks -query_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.cds.fa" -gff="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.gff" -hit_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potra.cds.fa" -spe=6 -blast2seq="bl2seq" -program="codeml" -blat="blat" -kaks_file="$OutDir/potra_potri.kaks.blat.txt" -detail="$OutDir/potra_potri.kaks.blat.detail.txt"


