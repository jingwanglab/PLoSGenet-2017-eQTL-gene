#! usr/bin/perl 

###perl vcf_maf.eqtls.pl /proj/b2011141/nobackup/eQTL_paper/vcf/freq2/maf/SwAsp_94_snps.freq_0_05.maf /proj/b2011141/nobackup/eQTL_paper/eQTLs/swasp_eqtls.txt

use strict;
use File::Basename;

my $frq_file=basename($ARGV[0],"\.maf");
my $dirname=dirname($ARGV[0]);
my $frq_annotation=$frq_file.".eqtls.maf";
my $Output=join "/", ($dirname,$frq_annotation);


open FREQ, "<", $ARGV[0] or die "cannot open the POS file: $!";
open EQTL, "<", $ARGV[1] or die "cannot open the POS file: $!";
open OUT, ">", $Output or die "cannot produce the OUT file: $!";

my %snp;

while(<FREQ>) {
	chomp;
        if ($_ =~ /SNP/) {next;}
	my @line=split(/\t/,$_);
	$snp{$line[0]}=join("\t",$line[1],$line[2]);
}

while (<EQTL>) {
	chomp;
	if ($_ =~ /eSNP/) {print OUT "$_\tMAF\tMAF_class\n";next;}
	my @line=split(/\t/,$_);
	print OUT "$_\t$snp{$line[0]}\n";
	
}
	
