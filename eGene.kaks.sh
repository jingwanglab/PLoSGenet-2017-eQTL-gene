#! /bin/bash -l

#SBATCH -A b2010014
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o eGene_gkaks.out
#SBATCH -e eGene_gkaks.err
#SBATCH -J eGene_gkaks.job
#SBATCH -t 2-00:00:00

module load bioinfo-tools
module load blat
module load blast
module load paml
##first only select the genes in the Potri gff3 file
##grep -v "gene" Ptrichocarpa_210_v3.0.gene.gff3 |cut -f 1 -d ";" |sed 's/.v3.0.*//g' |sed 's/ID=//g' > /proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.gff


#---------------Step 1: Using gKaKs to estimate the Ka/Ks for all the homologous genes defined in the study-----------------

OutDir="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/gkaks"
gkaks="/proj/b2011141/tools/gKaKs/gKaKs_v1.3.pl"

#perl $gkaks -query_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.cds.fa" -gff="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.gff" -hit_seq="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potra.cds.fa" -spe=6 -blast2seq="bl2seq" -program="codeml" -blat="blat" -sfbh=T -gene_pair="/proj/b2011141/nobackup/eQTL_paper/gene_list/Potra_Potri/Potri.Potra.txt" -kaks_file="$OutDir/potra_potri.kaks.txt" -detail="$OutDir/potra_potri.kaks.detail.txt" -min_length="50" -min_identity=0.7


#---------------Step 2: After modify the output file slightly, make the plot to separately compare the ka/ks for the four gene categories (core-eGene, core-noneGene, noncore-eGene, noncore-noneGenes)









