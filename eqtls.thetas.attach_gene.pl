#! usr/bin/perl 

#perl eqtls.thetas.attach_gene.pl /proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/eGene_core.thetas.txt  /proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/eGene_core.thetas.gz.gene.txt

use strict;
use File::Basename;

my $thetas=basename($ARGV[0],"\.thetas.txt");
my $dirname=dirname($ARGV[0]);
my $gene_thetas=$thetas.".thetas.gene.txt";
my $Output=join "/", ($dirname,$gene_thetas);


open THETAS, "<", $ARGV[0] or die "cannot open the pestPG file: $!";
open GENE, "<", $ARGV[1] or die "cannot open the Gene_list file: $!";
open OUT, ">", $Output or die "cannot produce the OUT file: $!";

my %gene_list;
while(<GENE>) {
	chomp;
	my @line=split(/\t/,$_);
	$gene_list{$line[0]}=$line[1];
}

while(<THETAS>) {
	chomp;
	if ($_ =~ /^Chr/) {print OUT "Gene\tChr\tPos\tnumSites\ttW.norm\ttP.norm\ttajD\tfulif\tfuliD\tfayH\tzengsE\n"; next;}
	my @line=split(/\t/,$_);
	my $chrom=join"_",($line[0],$line[1]);
	if (defined $gene_list{$chrom}){
		print OUT "$gene_list{$chrom}\t$_\n";
	}
}



