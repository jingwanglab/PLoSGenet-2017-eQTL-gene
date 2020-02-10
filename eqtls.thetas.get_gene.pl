#! usr/bin/perl 

#perl eqtls.thetas.get_gene.pl /proj/b2011141/nobackup/eQTL_paper/gene_list/thetas/eGene_core.thetas.gz.pestPG

use strict;
use File::Basename;

my $pestpg=basename($ARGV[0],"\.pestPG");
my $dirname=dirname($ARGV[0]);
my $gene_n=$pestpg.".gene.txt";
my $Output=join "/", ($dirname,$gene_n);


open PESTPG, "<", $ARGV[0] or die "cannot open the pestPG file: $!";
open OUT, ">", $Output or die "cannot produce the OUT file: $!";

my %gene_list;
my $gene;
while(<PESTPG>) {
	chomp;
        if ($_ =~ /#Potra/) {$gene = $_;next;}
        if ($_ =~ /Chr/) {next;}
	my @line=split(/\t/,$_);
        my $chrom=join"_",($line[1],$line[2]);
	$gene_list{$chrom}=$gene;

	print OUT "$chrom\t$gene_list{$chrom}\n";
}


