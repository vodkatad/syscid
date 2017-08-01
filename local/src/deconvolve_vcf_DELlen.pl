#!/usr/bin/env perl

use strict;
use warnings;
use lib $ENV{'BIOINFO_ROOT'}.'/prj/tba_variation/local/lib/Vcf';
use Vcf;
#use Data::Dump qw(dump);#sudo apt-get install libdata-dump-perl
use Getopt::Long;

$SIG{__WARN__} = sub {die @_};

my $usage=<<END;
Usage: $0 snp_info_file [-p] < vcf_file > sample_info_file
	-p|ignore_phase		continua anche se trova dei dati non fasati

L'ordine con cui compaiono le varianti nel sample_info_file e` consistente con la fase: 
controllato manualmente su un esempio, vedere local/share/docs/note_elena

snp_info_file contiene una sorta di indice degli SNP.
In stdout riporta, per ogni individuo e per ogni snp, le varianti sui due alleli (ordinati perche` assumiamo dati fasati).

.META: snp_info_file
	1	SNP_ID	snp_22_16057210
	2	ordered_SNP_ID	2
	3	chr	22
	4	pos	16057210
	5	ref	C
	6	alt	T

.META: stdout
	1	ordered_SNP_ID		2
	2	SAMPLE_ID	HG00105
	3	allele1		T
	4	allele2		C

Esempio:
	 /rogue_bis/data/bioinfotree/prj/tba_variation/local/bin/deconvolve_vcf -p snp_index > sample_variant </rogue_bis/data/bioinfotree/prj/tba_variation/local/share/data/piccolo.vcf
END
 
my $help=0;
my $ignore_phase=0;
GetOptions (
	'h|help' => \$help,
	'p|ignore_phase' => \$ignore_phase,
) or die($usage);

if($help){
	print $usage;
	exit(0);
}

my $snp_info_file = $ARGV[0];
if ($#ARGV != 0) {
	die($usage);
}

open SNP_INFO, ">$ARGV[0]";
my $vcf = Vcf->new(fh=>\*STDIN);
$vcf->parse_header();
my $ordered_id = 1;
while (my $x=$vcf->next_data_hash()) {
	my $chr = $x->{"CHROM"};
	my $coord = $x->{"POS"};
	my $reference = $x->{"REF"};
	my $alt = $x->{"ALT"};
	my $id = $x->{"ID"};
	print SNP_INFO "$id\t$ordered_id\t$chr\t$coord\t$reference";
	my $len = 0;
	foreach my $a ( @{ $alt } ) {
		if ($a eq "<DEL>") {	
			my $end = $x->{"INFO"}{"END"};
			if (defined($end)) {
				$len = $end-$coord;
				print SNP_INFO "\tDEL_$len";
			} else {
				die "I found a SNP with a DEL allele and no END! ", $id, "\n";	
				# Add an option to allow these cases? Is it reasonable?
				# It is possible to skip this SNP for these sample but the foreach has
				# to be changed.
			}
		} else {
			print SNP_INFO "\t$a";
		}
	}
	print SNP_INFO "\n";
	my @samples = keys %{ $$x{"gtypes"} }; 
	# sort is not needed here as long as long as we always print ordered_id and sample 
	# and then consider those to have sorted data.
	foreach my $sample (@samples) {
       		my ($alleles,$seps,$is_phased,$is_empty) = $vcf->parse_haplotype($x,$sample);
		die ("Phased data required.") if not $ignore_phase and not $is_phased;
		print "$ordered_id\t$sample";
		foreach my $a ( @{ $alleles } ) {
			if ($a eq "<DEL>") {	
				print "\tDEL_$len";
			} else {
				print "\t$a";
			}
		}
		print "\n";
	}	
	$ordered_id++;
}
close SNP_INFO;
