#! usr/bin/perl 

###perl vcf_maf.pl /proj/b2011141/nobackup/eQTL_paper/vcf/freq2/SwAsp_94_snps.freq_0_05.frq

use strict;
use File::Basename;

my $frq_file=basename($ARGV[0],"\.frq");
my $dirname=dirname($ARGV[0]);
my $frq_annotation=$frq_file.".maf";
my $Output=join "/", ($dirname,$frq_annotation);


open FREQ, "<", $ARGV[0] or die "cannot open the POS file: $!";
open OUT, ">", $Output or die "cannot produce the OUT file: $!";

while(<FREQ>) {
	chomp;
	my $maf;
        if ($_ =~ /CHROM/) {print OUT "SNP\tMAF\tMAF_class\n";next;}
	my @line=split(/\t/,$_);
        my $snp=join":",($line[0],$line[1]);
        if ($line[4]>=0.5) {
        $maf=$line[5];
	}
        else{$maf=$line[4];}
	
	my $maf_class="";
	if ($maf >= 0 && $maf <= 0.05) {
		$maf_class="0.00-0.05";
	}
	elsif ($maf > 0.05 && $maf <= 0.1) {
		$maf_class="0.05-0.10";
	}
	elsif ($maf > 0.1 && $maf <= 0.15) {
		$maf_class="0.10-0.15";
	}
	elsif ($maf > 0.15 && $maf <= 0.2) {
		$maf_class="0.15-0.20";
	}
	elsif ($maf > 0.2 && $maf <= 0.25) {
		$maf_class="0.20-0.25";
	}
	elsif ($maf > 0.25 && $maf <= 0.3) {
		$maf_class="0.25-0.30";
	}
	elsif ($maf > 0.3 && $maf <= 0.35) {
		$maf_class="0.30-0.35";
	}
	elsif ($maf > 0.35 && $maf <= 0.4) {
		$maf_class="0.35-0.40";
	}
	elsif ($maf > 0.4 && $maf <= 0.45) {
		$maf_class="0.40-0.45";
	}
	elsif ($maf > 0.45 && $maf <=0.5) {
		$maf_class="0.45-0.50";
	}

	print OUT "$snp\t$maf\t$maf_class\n";
}


